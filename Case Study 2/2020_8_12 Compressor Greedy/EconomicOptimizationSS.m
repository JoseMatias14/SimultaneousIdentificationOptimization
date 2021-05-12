function uOptk = EconomicOptimizationSS(mod,alg,x_var,u_var,f_var,xk_1,uk,fk,par)
% Apply MA optimization
%
% Inputs:
%    dxk = initial differential states
%    zk = initial algebraic states
%    uk = inputs
%    par = model parameters
%    yValuePlant = plant measurements
%    gradYPlantHat = plant gradient estimates
%    epsk,epsk = modifiers
%    M - inverse of the covariance matrix
%    H = output map function (H*zk = yk)
%    ma = method parameters
%
% Outputs:
%    uOptk = new setpoint
%    costk = objective function value
%    epsk,lambdak = updated modifiers

% Other m-files required: EconomicOptimizationSS.m, EconomicOptimizationSSCov.m, WellModel.m
% Subfunctions: OptimizationBoundsGasLiftRiser.m
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br
% January 2018; Last revision: 27-Feb-2018

import casadi.*

%[~,~,alg,x_var,u_var,f_var] = SystemModel(xk_1,uk,fk,par,modelFlag);

% Building NLP
% decision variables
w = {};
w = {w{:},x_var,u_var,f_var};

% constraints
g = {};
%Add the system model as constraints
g = {g{:},alg};

%Add operational constraints on P_out and q_n
g = {g{:},(x_var(5) + mod(5))*1e7};
g = {g{:},(x_var(3) + mod(3))/u_var(2)};

% objective function
%J = -((x_var(3) + mod(3)))/((x_var(8) + mod(6))*1e3);
J = -(100*(x_var(3) + mod(3)) + 100*(x_var(9) + mod(7)) - 100*(x_var(8) + mod(6)) - (x_var(10) + mod(8)));

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

% Assign solver
options = struct;
options.ipopt.print_level = 5;
%    options.ipopt.return_status = true;
F = nlpsol('solver','ipopt',nlp,options);

% Solving
%guess for the varibles
%setting lower and upper bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas(par);

%specifying values
wk = [];
lbw = [];
ubw = [];

lbw = [lbw;lbx;lbu;fk];%one is a dummy variable - represent time
ubw = [ubw;ubx;ubu;fk];
wk = [wk;xk_1;uk;fk];

lbg = [];
ubg = [];
lbg = [lbg;zeros(10,1);par.Pcomp_min;par.qN_min];%100 bar
ubg = [ubg;zeros(10,1);inf;par.qN_max];

% Solve
sol = F('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

fs = F.stats();
%clc

if fs.success

    % Extract Solution
    xSol = full(sol.x);
    %extracting the results (from symbolic to numerical)
    uOptk = xSol(11:12);

%     %apply filter to the optimal inputs
%     uk1 = uk + ma.Kinput*(uOptk - uk);

else
    disp('error');
    return
end
    
end
