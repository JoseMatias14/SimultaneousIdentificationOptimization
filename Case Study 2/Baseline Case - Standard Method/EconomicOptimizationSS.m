function [uOptk,costOFk,costOFkMod] = EconomicOptimizationSS(x0,u0,f0,par,alg,x_var,u_var,f_var,mod)
% Optimize the steady state model based on an economic criterion using MAy

% Inputs:
%    x0,u0,f0 - initial values for the states, inputs. Feed is fixed
%    par - system parameters
%    alg - symbolic model (algebraic equations)
%    x_var,u_var,f_var - symbolic states and inputs
%    mod - output modifier
%
% Outputs:
%    uOptk - optimal inputs
%    costOFk - predicted value of the OF using the non-modified model
%    costOFkMod - objective function (cost) - modified model
%
% Subfunctions: OptimizationBoundsSubseaGas.m

import casadi.*

%% Optimization 
%states and inputs bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas(par);

% decision variables
w = {};
wk = [];
lbw = [];
ubw = [];

% constraints
g = {};
lbg = [];
ubg = [];

% assigning the variables
w = {w{:},x_var,u_var,f_var};
% bounds
lbw = [lbw;lbx;lbu;f0];
ubw = [ubw;ubx;ubu;f0];
% initial guess
wk = [wk;x0;u0;f0];

% add the system model as constraints
g = {g{:},alg};
%Add operational constraints on P_out and q_n
g = {g{:},(x_var(5) + mod(5))*1e7};
g = {g{:},(x_var(3) + mod(3))/u_var(2)};
% bounds
lbg = [lbg;zeros(10,1);par.Pcomp_min;par.qN_min];
ubg = [ubg;zeros(10,1);inf;par.qN_max];
% Note: Equality constraint is defined by setting ub = lb with the same
% value

% declaring modified objective function
J = -((x_var(3) + mod(3)))/((x_var(8) + mod(6))*1e3);

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

% Assign solver
options = struct;
options.ipopt.print_level = 5;
F = nlpsol('solver','ipopt',nlp,options);

% Solve
sol = F('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

fs = F.stats();
%clc

if fs.success
    % Extract Solution
    xSol = full(sol.x);
    %optimal inputs
    uOptk = xSol(11:12);
    %non modified prediction of the OF
    costOFk = xSol(3)/(xSol(8)*1e3);
    %modified prediction of the OF
    costOFkMod = -full(sol.f);
    
else
    disp('error');
    return
end
    
end
