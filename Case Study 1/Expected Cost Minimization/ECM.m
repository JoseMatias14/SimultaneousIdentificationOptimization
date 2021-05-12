function [pik,ukArray,rho_k,pi_k_index,epsk,lambdak] = ECM(dxk_1,uk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,rho_k,H,ma)
% Apply ECM method (see Section 4.1 | Expected Cost Minimization)
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
%    pik = chosen control input
%    ukArray = new inputs setpoint (for all models in the set)
%    rho_k = updated belief vector
%    pi_k_index =  model with the highest belief
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, ReactorModel.m, OptimizationBoundsBlockReactor.m
% MAT-files required: covLow.mat
import casadi.*

%Measurement covariance - previously calculated
load('covLow');

%Computing inputs setpoint and belief for all models
ukArray = [];
probkArray = [];

for jj = 1:ma.nModels
    mod{jj} = [];
end

%updating modifiers and computing the optimal for the model
for jj = 1:ma.nModels
    
    %indicating which model should be used in this iteration
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %getting model predictions at current operational point
    [diff,x_var,u_var,L,~,measValue,gradMeas] = ReactorModel(dxk_1,uk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(H*yValuePlant - H*measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - H*gradMeas)';
    
    %computing the measurement modification (MAy)
    mod{jj} = H'*(epsk{jj} + lambdak{jj}'*(u_var - uk_1));
    
    %Steady state optimization of the modified problem
    [uOptk,~,~] = EconomicOptimizationSS(dxk_1,uk_1,par,diff,x_var,u_var,L,mod{jj});
    
    %store values
    ukArray = [ukArray; uOptk];

    % computed the measurement and gradient residuals and update the
    % probability of that the residuals como from a normal distribution
    res_k = [H*measValue - H*yValuePlant; H*gradMeas - gradYPlantHat];
    
    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    probkArray = [probkArray; mvnpdf(sqrt(1e1)*res_k,[],1e1*sigmaY)];
    
    clc

end

% normalizing cost array
probNArray = probkArray/(sum(probkArray));

%Bayesian update of the belief
rho_k = 1/(probNArray'*rho_k)*probNArray.*rho_k;

% find uk that minimizes cost function
VCost = zeros(ma.nModels,1);

for ii = 1:ma.nModels
    for jj = 1:ma.nModels

        %indicating which model should be used
        flag = zeros(1,ma.nModels);
        flag(jj) = 1;

        %reading model data evaluated at the uk's computed with all models
        [~,~,~,~,~,measValue,~] = ReactorModel(dxk_1,ukArray(ii),par,flag);

        %Modified bjective function (optimization: maximize C)
        mod_jj = H'*(epsk{jj} + lambdak{jj}'*(ukArray(ii) - uk_1));
        % we want to maximize: J := measValue(3) + [0, 0, 1, 0]*mod_jj
        % or minimize -J
        
        % building the expected cost expression
        VCost(ii,1) = VCost(ii,1) +  rho_k(jj)*(-measValue(3) - [0, 0, 1, 0]*mod_jj);

    end
end

    %choosing the policy that minimize the expected total cost
    [~,pi_k_index] = min(VCost);
    pik = ukArray(pi_k_index);


end

