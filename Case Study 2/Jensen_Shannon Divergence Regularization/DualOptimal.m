function [pik,rho_k,uk_lambda,epsk,lambdak] = DualOptimal(xk_1,uk_1,fk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,rho_k,H,ma,results)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create array of lambdas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ex: nint = 3
% Drawn one sample from each interval (0 - 1/3 | 1/3 - 2/3 | 2/3 - 1) based on an uniform distribution
l_lower = (0:1/par.nint:(1 - 1/par.nint))';
l_upper = (1/par.nint:1/par.nint:1)';

l_Array = l_lower + (l_upper-l_lower).*rand(par.nint,1);


if strcmp(results,'dual')
    % adding only the OF and only KL
    l_Array = [0; l_Array; 1];
elseif strcmp(results,'dual_005')
    l_Array = [0; 0.05; 1];
elseif strcmp(results,'dual_095')
    l_Array = [0; 0.95; 1];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. loop through all the models to build optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%states and inputs bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas(par);

% Building optimization problem
% Declaring variables
w = {};
wk = [];
lbw = [];
ubw = [];

% constraints
g = {};
lbg = [];
ubg = [];

%probability array 
probkArray = [];

%ma.nModels = 1;
%updating modifiers and computing the optimal for the model
for jj = 1:ma.nModels
    
    %indicating which model should be used
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %reading model data
    [alg{jj},x_var{jj},u_var{jj},f_var{jj},~,~,measValue,gradMeas,DxDu1{jj},gradU1{jj},DxDu2{jj},gradU2{jj}] = SystemModel(xk_1,uk_1,fk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(yValuePlant - measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - par.H*gradMeas)';
    
    %computing the measurement modification (MAy)
    mod{jj} = epsk{jj} + lambdak{jj}'*(u_var{jj} - uk_1);
    
    %Steady state optimization of the modified problem
    %uOptk = EconomicOptimizationSS(mod{jj},alg{jj},x_var{jj},u_var{jj},f_var{jj},xk_1,uk_1,fk_1,par);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2.1. use model equations as equality constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ===================================
    %     Optimization problem
    % ===================================
    % decision variables
    w = {w{:},x_var{jj},u_var{jj},f_var{jj},gradU1{jj},gradU2{jj}};
    lbw = [lbw;lbx;lbu;fk_1;-100*ones(2*10,1)];
    ubw = [ubw;ubx;ubu;fk_1;100*ones(2*10,1)];
    wk = [wk;xk_1;uk_1;fk_1;gradMeas(:,1);gradMeas(:,2)]; %iniitial guess
    
    
    % constraints (model)
    g = {g{:},alg{jj},DxDu1{jj},DxDu2{jj}};
    lbg = [lbg;zeros(10,1);zeros(2*10,1)]; %SS - dif and alg == 0
    ubg = [ubg;zeros(10,1);zeros(2*10,1)];
    
    %Add operational constraints on P_out and q_n
    g = {g{:},(x_var{jj}(5) + mod{jj}(5))*1e7};
    g = {g{:},(x_var{jj}(3) + mod{jj}(3))/u_var{jj}(2)};
    lbg = [lbg;par.Pcomp_min;par.qN_min];%100 bar
    ubg = [ubg;inf;par.qN_max];

    %%%%%%%%%%%%%
    % Update probability
    %%%%%%%%%%%%%
    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    % 2 inputs, so I need to stack the values
    res_k = [yValuePlant - measValue;gradYPlantHat(:,1) - par.H*gradMeas(:,1);gradYPlantHat(:,2) - par.H*gradMeas(:,2)];
    probkArray = [probkArray; mvnpdf(sqrt(5*1e1)*res_k,[],5*1e1*covMeas)];
    %probkArray = [probkArray; mvnpdf(sqrt(1e3)*res_k(1:8),[],1e3*covMeas(1:8,1:8))];
    %probkArray = [probkArray; mvnpdf(sqrt(1e3)*res_k([6, 14, 22]),[],1e3*covMeas([6, 14, 22],[6, 14, 22]))];
   
    clc

end

%forcing the input of all models to be the same
for jj = 2:ma.nModels
    g = {g{:},u_var{1} - u_var{jj}};
    lbg = [lbg;0;0]; % two inputs
    ubg = [ubg;0;0];
    
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

% Computing the optimal input and the Vk for all lambda values (to compute the expectation) 
uk_lambda = [];
Vk_lambda = [];

for ii = 1:length(l_Array)
    % OF
    J = 0;
    D = 0;
    for jj = 1:ma.nModels
    %    2.2. use modified model OF to generate optimization OF
        %J = J + rho_k(jj)*(-((x_var{jj}(3) + mod{jj}(3)))/((x_var{jj}(8) + mod{jj}(6))*1e3));
        J = J + rho_k(jj)*(-(100*(x_var{jj}(3) + mod{jj}(3)) + 100*(x_var{jj}(9) + mod{jj}(7)) - 100*(x_var{jj}(8) + mod{jj}(6)) - (x_var{jj}(10) + mod{jj}(8))));

        Dtemp = 0;
        for kk = 1:ma.nModels
                %Dtemp = Dtemp + rho_k(kk)*((covMeas(3,3) + (L{jj} - L{kk})^2)/covMeas(3,3) - 1/2);
                resKL = [H*x_var{jj};H*gradU1{jj};H*gradU2{jj}] - [H*x_var{kk};H*gradU1{kk};H*gradU2{kk}];
                Dtemp = Dtemp + rho_k(kk)*(1/2*(resKL'*inv(covMeas)*resKL));
        end
        D = D + rho_k(jj)*Dtemp;
    end
        
    % (Eq. 12)
    %OF = 0.7*J + 0.3*D;
    OF = l_Array(ii)*J + 1e-4*(1 - l_Array(ii))*(-D);
   
    %solve optimization problem
    % formalize it into an NLP problem
    nlp = struct('x',vertcat(w{:}),'f',OF,'g',vertcat(g{:}));
    
    % Assign solver
    options = struct;
    options.ipopt.print_level = 5;
    Vk = nlpsol('solver','ipopt',nlp,options);
    
    % ===================================
    %     Solving
    % ===================================
    % Solve
    sol = Vk('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    
    % catch error
    if Vk.stats.success ~=1
        msg = 'Error in the optimization';
        error(msg);
       
    else
        % if opt is found, update vectors
        uk_lambda = [uk_lambda, [full(sol.x(11));full(sol.x(12))]]; % you can use any u
        Vk_lambda = [Vk_lambda, full(sol.f)];
    end
    clc %clean ipopt output
    
end

    %choosing the policy that minimize the total cost (excluding the first
    %and last point (lambda = 0 and 1)
    [~,Vk_index] = min(Vk_lambda(2:end - 1));
    pik = uk_lambda(:,Vk_index + 1);


end

