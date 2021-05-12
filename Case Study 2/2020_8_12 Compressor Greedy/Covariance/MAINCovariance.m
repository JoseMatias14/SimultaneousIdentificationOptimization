% Model Structure Selection implementation in Subsea Compression System
% (gradient estimated via CDA - Central difference approximation)

% computes the covariance matrix

% Other m-files required:
%   InitialConditionSubseaGas.m
%   OptimizationBoundsCompressor.m
%   OptimizationBoundsSubseaGas.m
%   PlantModel.m

% MAT-files required:

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
results = 'covMatrix_high';

% For reproducibility
rng('default')  

%% Simulation tuning
%simulation length
tEnd = 300;

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
[xk,~,fk] = InitialConditionSubseaGas(par);

%(Forward) Finite Difference Approach - FDA
h = 0.05;%perturbation size

%% Simulation
count = 1;

% near the optimum
uk = [0.70; 0.74];

% Plant Array -- [y, yp(1), yp(2)]
yPlantArray = [];%plant measurements

while count <= tEnd

    fprintf('     iter. >>> %0.0f \n',count)

    %1. simulate plant (SS) and store plant data
    [xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);
    
    %2. Gradient estimation using CDA
    [gradPlantHat,uk_h1,uk_h2,yk_h1,yk_h2] = CDAGradient(H*xmk,xk,uk,fk,parPlant,h,H,count);

    yPlantArray = [yPlantArray,[H*xmk;gradPlantHat(:,1);gradPlantHat(:,2)]];

    %loop
    count = count + 1;

end

%% Computing variances
%variance theoretical 
covMeasT = diag([ones(8,1)*0.01;ones(16,1)*0.001].^2);

% computing covariance
covMeas = cov(yPlantArray');

covY = covMeas(1:8,1:8);
covGY = covMeas(8+1:end,8+1:end);

covYT = covMeasT(1:8,1:8);
covGYT = covMeasT(8+1:end,8+1:end);

save(results,'covMeas','covY','covGY','covMeasT','covYT','covGYT');