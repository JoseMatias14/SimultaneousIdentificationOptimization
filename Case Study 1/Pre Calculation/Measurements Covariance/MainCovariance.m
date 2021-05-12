% Computing the covariance to calculate the bayesian updated of the belief
% vector

% Other m-files required: 
% InitialConditionBlockReactor.m;
% OptimizationBoundsBlockReactor.m; ParametersBlockReactor.m; PlantModel;

% MAT-files required: none
clear  
close all
clc

%adding casadi
import casadi.*

% For reproducibility
rng('default')  

%saving file
results = 'covMatrix';

%%
nIter = 300; %number of experimental points for estimating covariance

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

%% Initialization
count = 1; %counter for the number of executions

% random inputs to compute the noise -- all inputs between 2.5 and 0.5
% uniform distribution
uk = (8 - 0.5).*rand + 0.5;

% Plant Array -- [y, yp(1), yp(2)]
yPlantArray = [];%plant measurements

%% Simulation
while count <= nIter
    
    fprintf('     iter. >>> %0.0f \n',count) 
    
    %1. simulate plant (SS) and store plant data
    [~,~,~,~,~,yValuePlant,FendGradMeas] = PlantModel(dx0,uk,parPlant);
    
    yPlantArray = [yPlantArray,[H*yValuePlant;H*FendGradMeas]];
     
    %Loop
    count = count + 1;
      
end

%% Computing variances
%variance theoretical 
covMeasT = diag([ones(3,1)*0.02;ones(3,1)*0.002].^2);
covYT = covMeasT(1:3,1:3);
covGYT = covMeasT(3+1:end,3+1:end);

% computing covariance from data
covMeas = cov(yPlantArray');
covY = covMeas(1:3,1:3);
covGY = covMeas(3+1:end,3+1:end);

save(results,'covMeas','covY','covGY','covMeasT','covYT','covGYT');