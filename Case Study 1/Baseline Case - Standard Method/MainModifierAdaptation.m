%Modifier Adaptation implementation in CSTR reactor  (gradient estimated via CDA - central difference approximation)
% Model structure estimation using Multiarmed bandit problem allocation strategies

% Three different models are available. Difference between the models lies in the reaction kinetic parameter &
% set of reactions (structural mismatch)
%
% Other m-files required:

% MAT-files required: PlantSurface_BR.mat (from:...\2018_03_13\PlantSurface)

% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.nno
% July 2019; Last revision: 29-Jul-2019

clear all 
close all
clc

%adding casadi
addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

% For reproducibility
rng('default')  

%saving file
results = 'results_standard';

%For plotting - previously calculated plant surface and optimum
load('PlantSurface_BR');

%%
nIter = 50;

%% Simulation tuning
%paramters
par = ParametersBlockReactor;
%initial condition
[dx0,u0] = InitialConditionBlockReactor;
[~,~,~,ubu] = OptimizationBoundsBlockReactor;

%System parameters for plant model
% (in case we want to add parametric uncertainty to the model, then par ~= parPlant)
parPlant = par;

%state output mapping (Xb, Xc & Xd)
H = [0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];
par.H = H;

%% MA tuning
%applying modifiers only to the cost function
ma.nMeas = 3;
ma.nInput = 1; %only Fa

%first order filters for the modifiers 
ma.Keps = 0.5*eye(ma.nMeas);%bias filter (same for all measurements) 0.9
ma.Klam = 0.5*eye(ma.nInput);%gradient filter 0.9
ma.Kinput = 0.25*eye(ma.nInput); %input filter 0.75

%% Multimodel tuning
%flag to choose between the different models available - structural plant
%model mismatch
ma.nModels = par.nr; %number of models
%Multi model structure
% m1, m2, m3: A + B -k1-> C & 2B -k2-> D
% m1: k1 = 0.75, k2 = 1.5 (TRUE MODEL)
% m2: k1 = 0.9,  k2 = 1.3 
% m3: k1 = 0.6,  k2 = 1.5 

% m4, m5, m6: A + B -k3-> 2C & C + 2B -k4-> D
% m4: k3 = 0.15, k4 = 0.05 
% m5: k3 = 0.2,  k4 = 0.03 
% m6: k3 = 0.1,  k4 = 0.01 

% m7, m8, m9: A + B -k1-> C & D as impurity
% m1: k1 = 0.1
% m2: k1 = 0.2
% m3: k1 = 0.3 

%% Gradient Estimation tuning
%(Forward) Finite Difference Approach - FDA
h = 0.75;%perturbation size

%% Initialization
count = 1; %counter for the number of executions

%1. Simulation
dxk = dx0;
uk = u0;
modelChoicek = 9;
rho_k = 1/par.nr*ones(par.nr,1);

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
probModelArray = rho_k;
totalModifierArray = [];

%model based optimization - which model is used in each iteration
modelArray = modelChoicek; %which model is used to optimize - includes probing points for plotting
modelArrayProbability = rho_k; %prob of the model to be used.

%% Simulation
while count <= nIter
    
    fprintf('     iter. >>> %0.0f \n',count) 
    
    %1. simulate plant (SS) and store plant data
    [~,~,~,~,dxk,yValuePlant,FendGradMeas] = PlantModel(dxk,uk,parPlant);
    
        %%%%%%%%%%%%%%%%% plant information
        inputPlantArray = [inputPlantArray, uk];
        yPlantArray = [yPlantArray, H*yValuePlant];
        yNoNoiseArray = [yNoNoiseArray, H*dxk];

        gradTemp = H*FendGradMeas;
        for i = 1:ma.nMeas
            gradMeasPlantArray{i} = [gradMeasPlantArray{i},gradTemp(i,:)'];
        end
        %%%%%%%%%%%%%%%%%%
        
    %2. estimating gradient - CDA
    [gradPlantHat,uk_h,yk_h,xk_h] = CDAGradient(dxk,uk,parPlant,h,H);
    
        %%%%%%%%%%%%%%%%% storing for ploting (inputs for plant excitation)
        inputPlantArray = [inputPlantArray,uk_h];
        yPlantArray = [yPlantArray,yk_h]; 
        yNoNoiseArray = [yNoNoiseArray, xk_h];

        for i = 1:ma.nMeas
            gradPlantHatArray{i} = [gradPlantHatArray{i}, gradPlantHat(i,:)'];%plant measurements (pressure) gradients
        end
        
         modelArray = [modelArray,4,4]; %indicating that this is a probing trial for gradient estimation
        %%%%%%%%%%%%%%%%% 

    %3. Optimizing 
    [ukArray,totalModifierk,rho_k,epsk,lambdak] = MAOptimization(dxk,uk,par,yValuePlant,gradPlantHat,epsk,lambdak,H,ma);    
    %optimal decision of the three models given the current modifiers
    uModelArray = [uModelArray, ukArray];
    %total modifier related to the three models
    totalModifierArray = [totalModifierArray, totalModifierk];
    
    %%%%%%%%%%%%%%%%%%%
    % Model choice strategy
    %%%%%%%%%%%%%%%%%%%
    %choosing the model with the highest probability
    [~,modelChoicek] = min(totalModifierk);
    
    %choosing the optimal input compute when model == modelChoice 
    uStar = uModelArray(modelChoicek,end);

    %% implementing input filter
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
    modelArray = [modelArray,modelChoicek];
    modelArrayProbability = [modelArrayProbability, rho_k];
    
    %loop 
    count = count + 1;
      
end

%calculating the results for the last point
%1. simulate plant (SS) and store plant data
[~,~,~,~,dxk,yValuePlant,FendGradMeas] = PlantModel(dxk,uk,parPlant);

%%%%%%%%%%%%%%%%% plant information
inputPlantArray = [inputPlantArray, uk];
yPlantArray = [yPlantArray, H*yValuePlant];
yNoNoiseArray = [yNoNoiseArray, H*dxk];

%% Saving the results and ploting
save(results,'uOptArray','epsArray','lambdaArray','gradPlantHatArray','gradMeasPlantArray','yPlantArray','yNoNoiseArray','inputPlantArray','uModelArray','totalModifierArray','modelArray','modelArrayProbability','parPlant','ma','nIter');

%plotting the results
ResultsPlot