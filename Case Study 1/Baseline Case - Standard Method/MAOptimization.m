function [ukArray,totalModifierArray,rho_k,epsk,lambdak] = MAOptimization(dxk_1,uk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,H,ma)
% Apply baseline method (see [8])
%
% Inputs:
%    dxk_1 = current value of the differential states (only used as initial guess to find the SS)
%    uk_1 = current value of the inputs
%    par = model parameters
%    yValuePlant = current plant measurements
%    gradYPlantHat = current plant gradient estimates
%    epsk_1,lambdak_1 = modifiers (past iteration)
%    rho_k = current belief vector
%    H = output map function (H*xk = yk)
%    ma = method parameters
%
% Outputs:
%    ukArray = new inputs setpoint (for all models in the set)
%    totalModifierArray = total modifier computed for all models
%    rho_k = updated belief vector
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, ReactorModel.m, OptimizationBoundsBlockReactor.m
% MAT-files required: covLow.mat
import casadi.*


addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%For covariance - previously calculated
load('covMatrix');

%Computing inputs setpoint and belief for all models
ukArray = [];
probkArray = [];

totalModifierArray = [];

%updating modifiers and computing the optimal for the model
for jj = 1:ma.nModels
    
    %indicating which model should be used
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %getting model predictions at current operational point
    [diff,x_var,u_var,L,~,measValue,gradMeas] = ReactorModel(dxk_1,uk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(H*yValuePlant - H*measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - H*gradMeas)';
    
    %computing the measurement modification (MAy)
    mod = H'*(epsk{jj} + lambdak{jj}'*(u_var - uk_1));
    
    %Steady state optimization of the modified problem
    [uOptk,~,~] = EconomicOptimizationSS(dxk_1,uk_1,par,diff,x_var,u_var,L,mod);
    
    %store values
    ukArray = [ukArray; uOptk];

    % total modifier
    res_eps = H*measValue - H*yValuePlant; 
    res_lam = H*gradMeas - gradYPlantHat;
    totalModifierArray = [totalModifierArray; sqrt(res_eps'*res_eps)^2 + sqrt(res_lam'*res_lam)^2];

    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    res_k = [H*measValue - H*yValuePlant; H*gradMeas - gradYPlantHat];
    probkArray = [probkArray; mvnpdf(sqrt(1e1)*res_k,[],1e1*covMeas)];
    
    clc

end

% normalizing probability array
rho_k = probkArray/(sum(probkArray));



end

