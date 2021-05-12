function [pik,ukArray,rho_k,pi_k_index,epsk,lambdak] = MAPGreedy(xk_1,uk_1,fk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,rho_k,H,ma)
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
%    costkArray = reward (for all models in the set)
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, ReactorModel.m
% Subfunctions: OptimizationBoundsBlockReactor.m
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br
% March 2018; Last revision: 29-Jul-2019

import casadi.*

%For covariance - previously calculated
%from: ...\Covariance
load('covMatrix');

%for computing inputs setpoint and cost for all models
ukArray = [];
probkArray = [];

for jj = 1:ma.nModels
    mod{jj} = [];
end

%updating modifiers and computing the optimal for the model
for jj = 1:ma.nModels
    
    %indicating which model should be used
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %reading model data
    [alg,x_var,u_var,f_var,~,~,measValue,gradMeas] = SystemModel(xk_1,uk_1,fk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(yValuePlant - measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - gradMeas)';
    
    %computing the measurement modification (MAy)
    mod{jj} = epsk{jj} + lambdak{jj}'*(u_var - uk_1);
    
    %Steady state optimization of the modified problem
    uOptk = EconomicOptimizationSS(mod{jj},alg,x_var,u_var,f_var,xk_1,uk_1,fk_1,par);
    
    %store values
    ukArray = [ukArray, uOptk];

    %reward: total modifier
    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    % 2 inputs, so I need to stack the values
    res_k = [yValuePlant - measValue;gradYPlantHat(:,1) - gradMeas(:,1);gradYPlantHat(:,2) - gradMeas(:,2)];
    probkArray = [probkArray; mvnpdf(sqrt(5*1e1)*res_k,[],5*1e1*covMeas)];
    %probkArray = [probkArray; mvnpdf(sqrt(1e3)*res_k(1:8),[],1e3*covMeas(1:8,1:8))];
    %probkArray = [probkArray; mvnpdf(sqrt(1e3)*res_k([6, 14, 22]),[],1e3*covMeas([6, 14, 22],[6, 14, 22]))];
   
    clc

end

if sum(probkArray) > 1e-7
    % normalizing cost array
    probNArray = probkArray/(sum(probkArray));

    %Bayesian update of the model prob.
    rho_k = 1/(probNArray'*rho_k)*probNArray.*rho_k;

    %hack to avoid probability going to zero
    if rho_k(1) < 0.001
        rho_k(1) = 0.001;
    end
    if rho_k(2) < 0.001
        rho_k(2) = 0.001;
    end
    if rho_k(3) < 0.001
        rho_k(3) = 0.001;
    end

    %renormalizing
    rho_k = rho_k/(sum(rho_k));
end %else - don't update

% find uk that minimizes cost function
VCost = zeros(ma.nModels,1);

for ii = 1:ma.nModels
    for jj = 1:ma.nModels

        %indicating which model should be used
        flag = zeros(1,ma.nModels);
        flag(jj) = 1;

        %reading model data
        [~,~,~,~,~,xVal,~,~] = SystemModel(xk_1,ukArray(:,ii),fk_1,par,flag);

        %Objective function (optimization: maximize C)
        mod_jj = epsk{jj} + lambdak{jj}'*(ukArray(:,ii) - uk_1);
        
        % or minimize -J
        %VCost(ii,1) = VCost(ii,1) + rho_k(jj)*(-((xVal(3) + mod_jj(3)))/((xVal(8) + mod_jj(6))*1e3));
        VCost(ii,1) = VCost(ii,1) + rho_k(jj)*(-(100*(xVal(3) + mod_jj(3)) + 100*(xVal(9) + mod_jj(7)) - 100*(xVal(8) + mod_jj(6)) - (xVal(10) + mod_jj(8))));

    end
end

    %choosing the policy that minimize the total cost
    [~,pi_k_index] = min(VCost);
    
    pik = ukArray(:,pi_k_index);


end

