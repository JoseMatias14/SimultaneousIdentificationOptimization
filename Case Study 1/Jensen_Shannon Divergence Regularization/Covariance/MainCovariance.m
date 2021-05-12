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

% computing covariance
covMeas = cov(yPlantArray');

covY = covMeas(1:3,1:3);
covGY = covMeas(3+1:end,3+1:end);

covYT = covMeasT(1:3,1:3);
covGYT = covMeasT(3+1:end,3+1:end);

save(results,'covMeas','covY','covGY','covMeasT','covYT','covGYT');