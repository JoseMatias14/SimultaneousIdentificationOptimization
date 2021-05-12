% Model Structure Selection implementation in Subsea Compression System
% (gradient estimated via CDA - Central difference approximation)
% ESCAPE paper

% Other m-files required:
%   InitialConditionSubseaGas.m
%   OptimizationBoundsCompressor.m
%   OptimizationBoundsSubseaGas.m
%   PlantModel.m

% MAT-files required:
%   PlantNoise.mat (from:...\Models\NoiseGenerator)
%   PlantSurface.mat (from:...\Models\SubSeaGasModel\Ongoing\SimplifiedModel)
% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.nno
% October 2019

clear
close all
clc

%adding casadi
%addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%saving file
results = 'testPlant';

% For reproducibility
rng('default')  

%Previously computed noise from standard normal distribution
load('PlantSurface');

%% Simulation tuning
%simulation length
tEnd = 200;

%paramters
par = ParametersSubseaGas;
%bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas;

%System parameters for plant model
% (in case we want to add parametric uncertainty to the model, then par ~= parPlant)
parPlant = par;

%state output mapping (Xb, Xc & Xd)
H = zeros(8,10);
H(1,1) = 1; %choke outlet pressure: P_choke_out*1e7 [Pa]
H(2,2) = 1; %volumetric flow from the reservoir: q_res [m3/s]
H(3,3) = 1; %inlet compressor volumetric flowrate: q_comp_in [m3/s]
H(4,4) = 1; %compressor outlet temperature: T_comp_out*1e2 [K]
H(5,5) = 1; %compressor outlet pressure: P_comp_out*1e7 [Pa]
H(6,8) = 1; %compressor power: Pow*1e4 [W]
H(7,9) = 1;%inlet pump volumetric flowrate: q_pump_in*1e-1 [m3/s]
H(8,10) = 1; %pump power: Pop*1e3 [W]
par.H = H;

%initial condition
[x0,u0,f0] = InitialConditionSubseaGas(par);

%% MSA tuning
%applying modifiers only to the cost function
ma.nMeas = 8;
ma.nInput = 2;

%first order filters for the modifiers
ma.Keps = 1*eye(ma.nMeas);%bias filter (same for all measurements)
ma.Klam = 1*eye(ma.nInput);%gradient filter
ma.Kinput = 0.2*eye(ma.nInput); %input filter

% Multimodel tuning
%flag to choose between the different models available - structural plant
%model mismatch
% if flag = [1,0,0] - model A (minRho = 1): A + B -> C & 2B -> D
% if flag = [0,1,0] - model B (minRho = 2): A + B -> C
% if flag = [0,0,1] - model C (minRho = 3): A + 2B -> C & B -> D
ma.modelType = [1,0,0;
                0,1,0;
                0,0,1];

ma.nModels = 3; %number of models

%% Gradient Estimation tuning
%(Forward) Finite Difference Approach - FDA
h = 0.05;%perturbation size

%% Initialization

%1. Simulation
xk = x0;
fk = f0;
uk = u0;%[u0, u0.*(0.9 + (1.05-0.9).*rand(2,9))];
modelChoicek = 1;
prob_k = 1/ma.nModels*ones(ma.nModels,1);


%2. MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsk{j} = zeros(ma.nMeas,1);
    lambdak{j} = zeros(ma.nInput,ma.nMeas);

end

%% For plotting - several values are stored for ploting
%1. Optimal values
uOptArray = u0;

%2. MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsArray{j} = epsk{j}; %dim OF = 1

    for i = 1:ma.nMeas
        lambdaArray{j,i} = lambdak{j}(:,i);
    end
end

%2. Gradient
for i = 1:ma.nMeas
    %plant measurements gradients
    gradMeasPlantArray{i} = [];
    %plant measurements gradients estimates
    gradPlantHatArray{i} = [];

end

%3. Plant
inputPlantArray = [];%optimal and points used for plant excitation and probing, i.e. all the inputs aplied during the simulation
yPlantArray = [];%plant measurements
yNoNoiseArray = [];%plant measurements no noise

%4. Models
% for computing the reward distribution
uModelArray = [];
probModelArray = prob_k;
uChosenModelArray = prob_k;

%model based optimization - which model is used in each iteration
modelArrayProb = modelChoicek; %which model has the highest probability
modelArrayU = modelChoicek;    %which model computed input is chosen

%% Simulation
count = 1;

while count <= tEnd

    fprintf('     iter. >>> %0.0f \n',count)

    %1. simulate plant (SS) and store plant data
    [xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);

        %%%%%%%%%%%%%%%%% plant information
        inputPlantArray = [inputPlantArray, uk];
        yPlantArray = [yPlantArray, H*xmk]; %H*xmk == yValuePlant

        gradTemp = H*gradPlant;
        for i = 1:ma.nMeas
            gradMeasPlantArray{i} = [gradMeasPlantArray{i},gradTemp(i,:)'];
        end
        %%%%%%%%%%%%%%%%%%

    %2. Gradient estimation using CDA
    [gradPlantHat,uk_h1,uk_h2,yk_h1,yk_h2] = CDAGradient(H*xmk,xk,uk,fk,parPlant,h,H,count);

        %%%%%%%%%%%%%%%%% plant information (inputs for plant excitation)
        inputPlantArray = [inputPlantArray,uk_h1,uk_h2];
        yPlantArray = [yPlantArray,yk_h1(1:ma.nMeas),yk_h2(1:ma.nMeas)];

        for i = 1:ma.nMeas
            gradPlantHatArray{i} = [gradPlantHatArray{i},gradPlantHat(i,:)'];
        end
        %%%%%%%%%%%%%%%%

   %3. Optimizing 
    [pi_k,ukArray,prob_k,pi_k_index,epsk,lambdak] = MAOptimization(xk,uk,fk,par,H*xmk,gradPlantHat,epsk,lambdak,prob_k,H,ma);    
    %optimal decision of the three models given the current modifiers
    uModelArray = [uModelArray, ukArray];
    %total modifier related to the three models
    probModelArray = [probModelArray, prob_k];
    % model-greedy controls that is implement at each period
    modelArrayU = [modelArrayU, pi_k_index];
    
    
    %%%%%%%%%%%%%%%%%%%
    % Model choice strategy
    %%%%%%%%%%%%%%%%%%%    
    %choosing the model with the highest probability
    [~,modelChoicek] = max(prob_k);
    
    %choosing the optimal input compute when model == modelChoice 
    uStar = pi_k;

    % implementing input filter
    uk = uk + ma.Kinput*(uStar - uk);
        
    %4. For plotting
    % Optimization
    uOptArray = [uOptArray, uk];
    
    % MA
    for j = 1:ma.nModels
        epsArray{j} = [epsArray{j}, epsk{j}];
        
        for i = 1:ma.nMeas
            lambdaArray{j,i} = [lambdaArray{j,i}, lambdak{j}(:,i)];
        end
    end
        
    % Model choice 
    modelArrayProb = [modelArrayProb,modelChoicek];
          
    %loop
    count = count + 1;

end

%calculating the results for the last point
%1. simulate plant (SS) and store plant data
[xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);

%%%%%%%%%%%%%%%%% plant information
inputPlantArray = [inputPlantArray, uk];
yPlantArray = [yPlantArray, H*xmk];

%% Saving the results and ploting
save(results,'inputPlantArray','yPlantArray','yNoNoiseArray','gradMeasPlantArray','gradPlantHatArray',...
              'uOptArray','uModelArray','probModelArray','uChosenModelArray','modelArrayProb','modelArrayU',...
              'ma','par','parPlant');

%% Plotting
inputs = {'u_{choke}','u_{comp}'};
measurements = {'Choke outlet pressure','Reservoir outflow','Compressor inflow','Compressor temperature','Compressor outlet pressure','Compressor efficiency','Compressor head','Compressor Power','Pump inflow','Pump power'};
units = {'P_{choke} [Pa]','q_{res} [m3/s]','q_{comp} [m3/s]','T_{comp} [K]','P_{comp} [Pa]','\nu [-]','H [m]','P_{ow} [W]','q_{pump} [m3/s]','P_{op} [W]'};
scaling = [1e7,1e0,1e0,1e2,1e7,1e0,1e1,1e4,1e-1,1e3];

%% 1. Input Comparison
leg = {'MAy','Opt'};

uSP = [];
for i = 1:count
    if i < 30
        uSP = [uSP, uStar_H];

    else
        if i < 60
            % degraded state 1
            uSP = [uSP, uStar_D1];

        else
            uSP = [uSP, uStar_D2];

        end
    end
end

f1 = figure(1);

    subplot(2,1,1,'FontSize',10)
        plot(1:count,uOptArray(1,:),'k',1:count,uSP(1,:),'k:')
        ylabel('u_{choke}','FontSize',10)
        xlabel('iteration','FontSize',10)
        xlim([1,count])
        ylim([0.7,1])
    subplot(2,1,2,'FontSize',10)
        plot(1:count,uOptArray(2,:),'b',1:count,uSP(2,:),'b:')
        ylabel('u_{comp}','FontSize',10)
        xlabel('iteration','FontSize',10) 
        xlim([1,count])
        ylim([0.7,1])

legend(leg,'Location','best','FontSize',9)
print(f1,'-r80','-dpdf','ResultsGradients.pdf');

%% 2. comparing the actual gradient value with the estimated
leg = {'plant','hat'};

f2 = figure(2);
for i = 1:ma.nMeas
    subplot(4,2,i,'FontSize',10)
        plot(1:count-1,gradMeasPlantArray{i}(1,:),'k:')
        hold on
        plot(1:count-1,gradPlantHatArray{i}(1,:),'k')

    ylabel(['g ',num2str(i)],'FontSize',10)

end
title('grad u choke')
legend(leg,'Location','best','FontSize',9)
print(f2,'-r80','-dpdf','ResultsGradients.pdf');

f3 = figure(3);
for i = 1:ma.nMeas
    subplot(4,2,i,'FontSize',10)
        plot(1:count-1,gradMeasPlantArray{i}(2,:),'b:')
        hold on
        plot(1:count-1,gradPlantHatArray{i}(2,:),'b')

    ylabel(['g ',num2str(i)],'FontSize',10)

end
title('grad u comp')
legend(leg,'Location','best','FontSize',9)
print(f3,'-r80','-dpdf','ResultsGradients.pdf');

%% 3. comparing choosen model with actual model
leg = {'plant','hat'};
mod = {'M_{1}','M_{2}','M_{3}'};

%real plant array
rpArray = [];
for count = 1:tEnd
    if count < 31 
        rpArray = [rpArray,1];
    elseif count < 91 && count > 60
        rpArray = [rpArray,2];
    elseif count > 121 
        rpArray = [rpArray,3];
    else %121 to 150
        rpArray = [rpArray,0];
    end
end

f4 = figure(4);
subplot(2,1,1)
    stairs(1:count,rpArray,'k:')
    hold on
    stairs(1:count,modelArrayProb(2:end),'k')

    ylabel('Model','FontSize',10)
    title('Model Choice')
    legend(leg,'Location','best','FontSize',9)


%plotting probability profile
subplot(2,1,2)
    area(0:count,probModelArray')
    
    xlim([0, count])
    ylim([0 1])  
    
    title('Product concentration (C)')
    xlabel('SS periods [-]')
    ylabel('Probability of model')
    legend(mod,'FontSize',6,'Location','best')

print(f4,'-r80','-dpdf','ResultsModels.pdf');
hold off


beep

% %% 4. comparing estimated parameters
% 
% %real plant array
% rpThetaArray{1} = [];
% rpThetaArray{2} = [];
% rpThetaArray{3} = [];
% 
% for count = 1:tEnd
%     if count < 31 
%         rpThetaArray{1} = [rpThetaArray{1},[0.582;-2.398;2.75;-3.969;4.303]];
%     elseif count < 91 && count > 60
%         rpThetaArray{2} = [rpThetaArray{2},[-0.5;1.8;-0.85]];
%     elseif count > 121
%         rpThetaArray{3} = [rpThetaArray{3},[0.43;0.7;5;1;2]];
%     else 
%         rpThetaArray{1} = [rpThetaArray{1},[0;0;0;0;0]];
%         rpThetaArray{2} = [rpThetaArray{2},[0;0;0]];
%         rpThetaArray{3} = [rpThetaArray{3},[0;0;0;0;0]];
%     end
% end
% 
% f5 = figure(5);
% color = ['k','b','r'];
% for i = 1:ma.nModels
%         plot(1:count,rpThetaArray{i},color(i),'LineSpec',':')
%         hold on
%         plot(1:count,thetaHat{i},color(i))
% 
% end
% 
% ylabel('Parameters','FontSize',10)
% title('Estimated parameters')
% print(f4,'-r80','-dpdf','ResultsParameters.pdf');
% 
% 
% 
