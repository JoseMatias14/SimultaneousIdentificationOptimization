function [pik,rho_k,uk_lambda,epsk,lambdak] = SystemOptimization(dxk_1,uk_1,par,yValuePlant,gradYPlantHat,epsk_1,lambdak_1,rho_k,H,ma)
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
% June 2020; Last revision: 

%%
%%%%%%%%%%%%%%%%%%%%%%
% CHANGE PATH HERE
%%%%%%%%%%%%%%%%%%%%%%
%addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%%

load('covMatrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create array of lambdas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ex: nint = 3
% Drawn one sample from each interval (0 - 1/3 | 1/3 - 2/3 | 2/3 - 1) based on an uniform distribution
l_lower = (0:1/par.nint:(1 - 1/par.nint))';
l_upper = (1/par.nint:1/par.nint:1)';

l_Array = l_lower + (l_upper-l_lower).*rand(par.nint,1);

% adding only the OF and only KL
l_Array = [0; l_Array; 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. loop through all the models to build optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%states and inputs bounds
[lbx,lbu,ubx,ubu] = OptimizationBoundsBlockReactor;

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

for jj = 1:ma.nModels
    
    %indicating which model should be used
    %Multi model structure
    % flag(1,1) = 1 -> m1: k1 = 0.75, k2 = 1.2
    % flag(1,2) = 1 -> m2: k1 = 0.72,  k2 = 1.5
    % flag(1,3) = 1 -> m3: k1 = 0.75,  k2 = 1.5  (TRUE MODEL)
    % flag(1,4) = 1 -> m4: k3 = 0.15, k4 = 0.05
    % flag(1,5) = 1 -> m5: k3 = 0.2,  k4 = 0.03
    % flag(1,6) = 1 -> m6: k3 = 0.1,  k4 = 0.01
    % flag(1,7) = 1 -> m1: k1 = 0.1
    % flag(1,8) = 1 -> m2: k1 = 0.2
    % flag(1,9) = 1 -> m3: k1 = 0.3
    flag = zeros(1,ma.nModels);
    flag(jj) = 1;
    
    %reading model data
    % Outputs:
    %   symbolic (CASADI): 
    %        diff := dC/dt
    %        x    := C
    %        u    := Fa
    %        L    := OF = -C(3) (minimize)
    %        DxDu := sensitivity equations
    %       gradU := dC/du
    %   numerical:
    %       measValue := C
    %       gradMeas  := dC/du
    [diff{jj},x_var{jj},u_var{jj},L{jj},DxDu{jj},gradU{jj},~,measValue,gradMeas] = ReactorModel(dxk_1,uk_1,par,flag);
    
    %updating the modifers with the model information at uk - filtered modifiers
    epsk{jj} = (eye(ma.nMeas) - ma.Keps)*epsk_1{jj} +  ma.Keps*(H*yValuePlant - H*measValue);
    lambdak{jj} = (eye(ma.nInput) - ma.Klam)*lambdak_1{jj} + ma.Klam*(gradYPlantHat - H*gradMeas)';
    
    %computing the measurement modification (MAy)
    mod{jj} = H'*(epsk{jj} + lambdak{jj}'*(u_var{jj} - uk_1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2.1. use model equations as equality constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ===================================
    %     Optimization problem
    % ===================================
    % decision variables
    w = {w{:},x_var{jj},u_var{jj},gradU{jj}};
    lbw = [lbw;lbx;lbu;-100*ones(length(lbx),1)];
    ubw = [ubw;ubx;ubu;100*ones(length(lbx),1)];
    wk = [wk;dxk_1;uk_1;gradMeas]; %iniitial guess

    % constraints (model)
    g = {g{:},diff{jj},DxDu{jj}};
    lbg = [lbg;zeros((length(dxk_1)),1);zeros((length(gradMeas)),1)]; %SS - dif and alg == 0
    ubg = [ubg;zeros((length(dxk_1)),1);zeros((length(gradMeas)),1)];
    
    % Add operational constraints (no operational constraints in this case)
    
    %reward: total modifier
    % returns the pdf of the normal distribution with mean mu and standard deviation sigma, evaluated at the values in x.
    res_k = [H*measValue - H*yValuePlant; H*gradMeas - gradYPlantHat];
    
    %%%%%%%%%%%%%
    % CHECK HERE
    %%%%%%%%%%%%%
    % I compute covMeas offline  (...\Covariance\MainCovariance.m)
    % I ran 300 times the "plant" using the same input. Then, compute the
    % covariance of the measurements (Cmeas and dCmeas/du)
    % Then, I calculated the probability of res_k coming out of a normal
    % distribution with mean zero [] and covariance covMeas. 
    % I scaled the values by 1e1 [Var(cX) = c^2*Var(X)]
    probkArray = [probkArray; mvnpdf(sqrt(1e1)*res_k,[],1e1*covMeas)];

end

%forcing the input of all models to be the same
for jj = 2:ma.nModels
    g = {g{:},u_var{1} - u_var{jj}};
    lbg = [lbg;0];
    ubg = [ubg;0];
    
end

%    2.3. compute pi (mean and cov-var matrix)
    %%%%%%%%%%%%%
    % ADD CODE HERE!!!
    %%%%%%%%%%%%%
    
% 3. with pi's compute JS
% normalizing cost array
probNArray = probkArray/(sum(probkArray));

%Bayesian update of the model prob. (Eq. 9)
rho_k = 1/(probNArray'*rho_k)*probNArray.*rho_k;

% Computing the optimal input and the Vk for all lambda values (to compute the expectation) 
uk_lambda = [];
Vk_lambda = [];

for ii = 1:length(l_Array)
    % OF
    J = 0;
    D = 0;
    for jj = 1:ma.nModels
    %    2.2. use modified model OF to generate optimization OF
        J = J + rho_k(jj)*(L{jj} + [0, 0, 1, 0]*mod{jj});
        Dtemp = 0;
        for kk = 1:ma.nModels
                %Dtemp = Dtemp + rho_k(kk)*((covMeas(3,3) + (L{jj} - L{kk})^2)/covMeas(3,3) - 1/2);
                resKL = [H*x_var{jj};H*gradU{jj}] - [H*x_var{kk};H*gradU{kk}];
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
        uk_lambda = [uk_lambda; full(sol.x(5))]; % you can use any u
        Vk_lambda = [Vk_lambda; -full(sol.f)];
    end
    clc %clean ipopt output
    
end

    %choosing the policy that minimize the total cost (excluding the first
    %and last point (lambda = 0 and 1)
    [~,Vk_index] = min(Vk_lambda(2:end - 1));
    pik = uk_lambda(Vk_index + 1);
    
end

