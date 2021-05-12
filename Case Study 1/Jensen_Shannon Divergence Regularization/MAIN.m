% Steady-state optimization via Modifier Adaptation of a CSTR reactor  
% (gradient estimated via CDA - central difference approximation)

% Nine different models are available. Difference between the models lies in the reaction kinetic parameter &
% set of reactions (structural and parametric mismatch)
%
% Other mat-files required:
% MAT-files required: PlantSurface_BR.mat (from:...\2018_03_13\PlantSurface)

% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.nno
% July 2019; Last revision: 05-Jun-2020

clear 
close all
clc

% For reproducibility
rng('default')  

%saving file
results = 'results_DUAL';

%% Simulation tuning
% loading model parameters
par = ParametersBlockReactor;

% loading initial condition
[dx0,u0] = InitialConditionBlockReactor;

% loading input bounds
[~,~,~,ubu] = OptimizationBoundsBlockReactor;

% System parameters for plant model
% (in case we want to add parametric uncertainty to the model, then par ~= parPlant)
parPlant = par;

%state output mapping (Xb, Xc & Xd)
H = [0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];
par.H = H;

% simulation d
nIter = 50;

%% MA tuning
%number of system measurements
ma.nMeas = 3;
%number of inputs
ma.nInput = 1; %only Fa

% filters for the modifiers 
ma.Keps = 0.5*eye(ma.nMeas); %zero-order modifier filter
ma.Klam = 0.5*eye(ma.nInput);%first-orde modifier filter
ma.Kinput = 0.25*eye(ma.nInput); %input filter

%% Multimodel tuning
%flag to choose between the different models available - structural plant
%model mismatch
ma.nModels = par.nr; %number of models
%Multi model structure
% m1, m2, m3: A + B -k1-> C & 2B -k2-> D
% m1: k1 = 0.75, k2 = 1.2 
% m2: k1 = 0.72,  k2 = 1.5 
% m3: k1 = 0.75,  k2 = 1.5  (TRUE MODEL)

% m4, m5, m6: A + B -k3-> 2C & C + 2B -k4-> D
% m4: k3 = 0.15, k4 = 0.05 
% m5: k3 = 0.2,  k4 = 0.03 
% m6: k3 = 0.1,  k4 = 0.01 

% m7, m8, m9: A + B -k1-> C & D as impurity
% m7: k1 = 0.1
% m8: k1 = 0.2
% m9: k1 = 0.3 

%% Gradient Estimation tuning
% Central Difference Approach - CDA
%perturbation size
h = 0.75;

%% Initialization
count = 1; %counter for the number of executions

%1. Simulation
dxk = dx0;
uk = u0;
modelChoicek = 9;

% uniform probability
prob_k = 1/par.nr*ones(par.nr,1);

%2. MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsk{j} = zeros(ma.nMeas,1); 
    lambdak{j} = zeros(ma.nInput,ma.nMeas); 

end

%% For plotting
%1. Optimal values
uOptArray = u0;

%2. MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsArray{j} = epsk{j}; 
    
    for i = 1:ma.nMeas
        lambdaArray{j,i} = lambdak{j}(:,i);
    end  
end

%2. Gradient 
for i = 1:ma.nMeas
    %plant gradients
    gradMeasPlantArray{i} = [];
    
    %plant gradients estimates
    gradPlantHatArray{i} = [];

end

%3. Plant
inputPlantArray = [];%optimal and points used for plant excitation and probing, i.e. all the inputs aplied during the simulation
yPlantArray = [];%plant measurements
yNoNoiseArray = [];%plant measurements no noise

%4. Models
% for computing the reward distribution
probModelArray = prob_k;

%model based optimization - which model is used in each iteration
modelArrayProb = modelChoicek; %which model has the highest probability

%store all the inputs used for finding the expected value
ukLambdaArray = [];%optimal and points used for plant excitation and probing, i.e. all the inputs aplied during the simulation


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
    
        %%%%%%%%%%%%%%%%% for plotting (plant info - plant excitation for gradient estimation)
        inputPlantArray = [inputPlantArray,uk_h];
        yPlantArray = [yPlantArray,yk_h]; 
        yNoNoiseArray = [yNoNoiseArray, xk_h];

        for i = 1:ma.nMeas
            gradPlantHatArray{i} = [gradPlantHatArray{i}, gradPlantHat(i,:)'];%plant measurements (pressure) gradients
        end
        
         modelArrayProb = [modelArrayProb,10,10]; %10 indicates that this is a probing trial for gradient estimation
        %%%%%%%%%%%%%%%%% 

    %3. Optimizing 
    [pi_k,prob_k,uk_lambda,epsk,lambdak] = SystemOptimization(dxk,uk,par,yValuePlant,gradPlantHat,epsk,lambdak,prob_k,H,ma);    

    
    % implementing input filter
    uk_lambda = uk + ma.Kinput*(uk_lambda - uk); %array for checking
    uk = uk + ma.Kinput*(pi_k - uk); % actual implemented value
    
        
    	%%%%%%%%%%%%%%%%% For plotting
        %total modifier related to the three models
        probModelArray = [probModelArray, prob_k];
        
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
        [~,modelChoicek] = max(prob_k);
        modelArrayProb = [modelArrayProb,modelChoicek];
        ukLambdaArray = [ukLambdaArray, uk_lambda];
        %%%%%%%%%%%%%%%%%
    
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
save(results,'uOptArray','epsArray','lambdaArray','gradPlantHatArray','gradMeasPlantArray','yPlantArray','yNoNoiseArray','inputPlantArray','ukLambdaArray','probModelArray','modelArrayProb','nIter','par','ma','parPlant','nIter');

%plotting the results
ResultsPlot