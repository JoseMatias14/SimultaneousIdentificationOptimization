function [ukArray,totalModifierArray,rho_k,epsk,lambdak] = MAOptimization(dxk_1,uk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,H,ma)
% Apply MA optimization
%
% Inputs:
%    dxk_1 = current value of the differential states (only used as initial guess to find the SS)
%    uk_1 = current value of the inputs
%    par = model parameters
%    yValuePlant = plant measurements
%    gradYPlantHat = plant gradient estimates
%    epsk_1,lambdak_1 = modifiers (past iteration)
%    H = output map function (H*xk = yk)
%    ma = method parameters
%
% Outputs:
%    ukArray = new inputs setpoint (for all models in the set)
%    totalModifierArray = totalModifier
%    rho_k = probability of the model being the correct one
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, ReactorModel.m
% Subfunctions: OptimizationBoundsBlockReactor.m
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br
% March 2018; Last revision: 29-Jul-2019

addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%For covariance - previously calculated
%from: \\home.ansatt.ntnu.no\joseoa\Documents\JOAM\Models\MultiArmedBanditApproach\2020_02_18\Low Noise
load('covLow');

%for computing inputs setpoint and cost for all models
ukArray = [];
totalModifierArray = [];
probkArray = [];

%updating modifiers and computing the optimal for the model
for j = 1:ma.nModels
    
    %indicating which model should be used
    flag = zeros(1,9);
    flag(j) = 1;
    
    %reading model data
    [diff,x_var,u_var,L,~,measValue,gradMeas] = ReactorModel(dxk_1,uk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{j} = (eye(ma.nMeas) - ma.Keps)*epsk_1{j} +  ma.Keps*(H*yValuePlant - H*measValue);
    lambdak{j} = (eye(ma.nInput) - ma.Klam)*lambdak_1{j} + ma.Klam*(gradYPlantHat - H*gradMeas)';
    
    %computing the measurement modification (MAy)
    mod = H'*(epsk{j} + lambdak{j}'*(u_var - uk_1));
    
    %Steady state optimization of the modified problem
    [uOptk,~,~] = EconomicOptimizationSS(dxk_1,uk_1,par,diff,x_var,u_var,L,mod);
    
    %store values
    ukArray = [ukArray; uOptk];

    %reward: total modifier
    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    res_eps = H*measValue - H*yValuePlant; 
    res_lam = H*gradMeas - gradYPlantHat;
    totalModifierArray = [totalModifierArray; sqrt(res_eps'*res_eps)^2 + sqrt(res_lam'*res_lam)^2];
    
    res_k = [H*measValue - H*yValuePlant; H*gradMeas - gradYPlantHat];
    probkArray = [probkArray; mvnpdf(sqrt(1e1)*res_k,[],1e1*sigmaY)];
    
    clc

end

% normalizing cost array
rho_k = probkArray/(sum(probkArray));



end

