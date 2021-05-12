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
results = 'baseline';

% For reproducibility
rng('default')  

%% Simulation tuning
%simulation length
tEnd = 250;

%paramters
par = ParametersSubseaGas;
%bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas;

%System parameters for plant model
% (in case we want to add parametric uncertainty to the model, then par ~= parPlant)
parPlant = par;

%state output mapping
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
ma.Keps = 0.75*eye(ma.nMeas);%bias filter (same for all measurements) 1
ma.Klam = 0.5*eye(ma.nInput);%gradient filter 1
ma.Kinput = 0.25*eye(ma.nInput); %input filter

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
%generating square wave - auxiliary
tAux = 1:tEnd;
sqAux = square(0.1*tAux);

%2. MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsk{j} = zeros(ma.nMeas,1);
    lambdak{j} = zeros(ma.nInput,ma.nMeas);

end

% maintenance stop flags
% if we detect degraded 2 state for par.maintPeriods periods, we stop the
% process for 15 SS periods.
maintFlag = 0;
% indicating that the process is stopt 
stopFlag = 0;

% this correction term corrects the simulation counter to indicate what is
% happening on the plant. If no maintenance is performed, count =
% countPlant. However, if the compressor is repaired, the countPlant is set
% to zero to indicate a new compressor
counterCorrection = 0; 

% simulation counter
count = 1;
countPlant = count; % indicates the true plant state. It is re.initialized everytime the compressor is repaired

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

%model based optimization - which plant model is used in each iteration
plantModelArray = [];

% array checking if maintenance has been performed
maintenanceArray = [];

%% Simulation
while count <= tEnd

    fprintf('     iter. >>> %0.0f \n',count)
    %1. simulate plant (SS) and store plant data
    if sqAux(count) == 1
        fk(3) = 1.1;
    else %sqAux == -1
        fk(3) = 1.2;
    end
    
    if stopFlag == 0
        
        % changing the input
        [xk,xmk,gradPlant,pModel] = PlantModel(xk,uk,fk,parPlant,countPlant);
        
            %%%%%%%%%%%%%%%%% plant information
            inputPlantArray = [inputPlantArray, uk];
            yPlantArray = [yPlantArray, H*xmk]; %H*xmk == yValuePlant

            gradTemp = H*gradPlant;
            for i = 1:ma.nMeas
                gradMeasPlantArray{i} = [gradMeasPlantArray{i},gradTemp(i,:)'];
            end
            %%%%%%%%%%%%%%%%%%
        
        %2. Gradient estimation using CDA
        [gradPlantHat,uk_h1,uk_h2,yk_h1,yk_h2] = CDAGradient(H*xmk,xk,uk,fk,parPlant,h,H,countPlant);
        
            %%%%%%%%%%%%%%%%% plant information (inputs for plant excitation)
            inputPlantArray = [inputPlantArray,uk_h1,uk_h2];
            yPlantArray = [yPlantArray,yk_h1(1:ma.nMeas),yk_h2(1:ma.nMeas)];

            for i = 1:ma.nMeas
                gradPlantHatArray{i} = [gradPlantHatArray{i},gradPlantHat(i,:)'];
            end
            %%%%%%%%%%%%%%%%
        
        %3. Optimizing
        [pi_k,~,~,pi_k_index,epsk,lambdak] = OnlineMaintenance(xk,uk,fk,par,H*xmk,gradPlantHat,epsk,lambdak,prob_k,H,ma);
        
        %%%%%%%%%%%%%%%%%%%
        % Model choice strategy
        %%%%%%%%%%%%%%%%%%%
        %choosing the model with the highest probability
        modelChoicek = pi_k_index;
        
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
        
        % Plant true model
        plantModelArray = [plantModelArray, pModel]; %maintenance
        
     else
        % this part of the code represents a maintenance stop
        
            %%%%%%%%%%%%%%%%% plant information
            inputPlantArray = [inputPlantArray, zeros(2,1),zeros(2,1),zeros(2,1)]; %optimization + probing
            yPlantArray = [yPlantArray, zeros(8,1), zeros(8,1), zeros(8,1)]; 

            for i = 1:ma.nMeas
                gradMeasPlantArray{i} = [gradMeasPlantArray{i}, zeros(1,2)'];
                gradPlantHatArray{i} = [gradPlantHatArray{i},zeros(1,2)'];
            end
            %%%%%%%%%%%%%%%%

        %4. For plotting
        % Optimization
        uOptArray = [uOptArray, zeros(2,1)];

        % MA
        for j = 1:ma.nModels
            epsArray{j} = [epsArray{j}, zeros(8,1)];

            for i = 1:ma.nMeas
                lambdaArray{j,i} = [lambdaArray{j,i}, zeros(2,1)];
            end
        end

        % Model choice 
        modelArrayProb = [modelArrayProb,-1];
        
        % Plant true model
        plantModelArray = [plantModelArray, 0]; %maintenance
        
    end
    
     % updating maintenance flag
    if modelChoicek == 3
        maintFlag = maintFlag + 1;
    end
    
    % maintanence stop
    if maintFlag > par.maintPeriods && maintFlag < par.stopPeriods
        if pModel == 3 % compressor is degraded
            stopFlag  = 1;
        else
            stopFlag  = 0;
        end
        
    elseif maintFlag >= par.stopPeriods
        stopFlag  = 0;
        maintFlag = 0;
        % indicating a new compressor
        counterCorrection = count;
    end
    
    %loop
    count = count + 1;
    
    if stopFlag == 0
        countPlant = count - counterCorrection;
    end
    
    % checking if maintenace has been performed
    if stopFlag == 0 && maintFlag <= par.maintPeriods
        maintValue = 0; %NO 
    elseif stopFlag == 0 && maintFlag > par.maintPeriods
        maintValue = 1; %ONLY CHECKING, NO MAINTENANCE 
        maintFlag = 0;
    else % stop flag == 1
        maintValue = 2; %YES 
    end
    
    maintenanceArray = [maintenanceArray, maintValue];
    
end

%calculating the results for the last point
%1. simulate plant (SS) and store plant data
[xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);

%%%%%%%%%%%%%%%%% plant information
inputPlantArray = [inputPlantArray, uk];
yPlantArray = [yPlantArray, H*xmk];

%% Saving the results and ploting
save(results,'inputPlantArray','yPlantArray','yNoNoiseArray','gradMeasPlantArray','gradPlantHatArray',...
              'uOptArray','uModelArray','probModelArray','uChosenModelArray','modelArrayProb','modelArrayU','plantModelArray','maintenanceArray',...
              'tEnd','ma','par','parPlant');

ResultsPlotting          
          
beep