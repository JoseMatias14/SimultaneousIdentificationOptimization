function [pik,ukArray,rho_k,pi_k_index,epsk,lambdak] = OnlineMaintenance(xk_1,uk_1,fk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,rho_k,H,ma)
% Apply baseline method (see [8])
%
% Inputs:
%    xk_1 = current value of the states (only used as initial guess to find the SS)
%    uk_1, fk_1 = current value of the inputs and feed
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
%    pik = inputs for the models with the smallest total modifier
%    rho_k = belief vector (not updated here)
%    pi_k_index = model with the smallest total modifier
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, SystemModel.m
% MAT-files required: none

import casadi.*

%for computing inputs setpoint and cost for all models
ukArray = [];
modArray = [];

for jj = 1:ma.nModels
    mod{jj} = [];
end

%updating modifiers and computing the optimal for the model
for jj = 1:ma.nModels
    
    %indicating which model should be used
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %getting model predictions at current operational point
    [alg,x_var,u_var,f_var,~,~,measValue,gradMeas] = SystemModel(xk_1,uk_1,fk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(yValuePlant - measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - gradMeas)';
    
    %computing the measurement modification (MAy)
    mod{jj} = epsk{jj} + lambdak{jj}'*(u_var - uk_1);
    
    %computing total modifier
    res1 = yValuePlant - measValue;
    res2 = (gradYPlantHat - gradMeas)';
    tMod = 0.5*(res1'* res1) + 0.5*(([res2(1,:),res2(2,:)])*([res2(1,:),res2(2,:)])');

    modArray = [modArray, tMod];
    
    %Steady state optimization of the modified problem
    [uOptk,~,~] = EconomicOptimizationSS(xk_1,uk_1,fk_1,par,alg,x_var,u_var,f_var,mod{jj});
    clc
    
    %store values
    ukArray = [ukArray, uOptk];

end
    %choosing the input (policy) that minimizes the total cost
    [~,pi_k_index] = min(modArray);
    
    pik = ukArray(:,pi_k_index);

end

