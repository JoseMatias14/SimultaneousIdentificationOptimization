function [uOptk,costOFk,costOFkMod] = EconomicOptimizationSS(dx0,u0,par,diff,x_var,u_var,L,mod)
% Optimize the steady state model based on an economic criterion using MAy

% Inputs:
%    dx0,u0 - initial guess for the states and inputs
%    par - system parameters
%    diff - symbolic model (differential equations)
%    x_var,u_var - symbolic states and inputs
%    L - cost function
%    mod - output modifier
%
% Outputs:
%    uOptk - optimal inputs
%    costOFk - predicted value of the OF using the non-modified model
%    costOFkMod - objective function (cost) - modified model
%
% Subfunctions: OptimizationBoundsBlockReactor.m

import casadi.*

%% Optimization 
%states and inputs bounds
[lbx,lbu,ubx,ubu] = OptimizationBoundsBlockReactor;

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
w = {w{:},x_var,u_var};
% bounds
lbw = [lbw;-inf;-inf;-inf;-inf;lbu];
ubw = [ubw;inf;inf;inf;inf;ubu];
% initial guess
wk = [wk;dx0;u0];

% add the system model as constraints
g = {g{:},diff};
lbg = [lbg;zeros(par.nc,1)]; %SS -> dif and alg == 0
ubg = [ubg;zeros(par.nc,1)];
% Note: Equality constraint is defined by setting ub = lb with the same
% value

% bounds on the modified states
g = {g{:},x_var + mod};
lbg = [lbg;lbx]; 
ubg = [ubg;ubx];

% declaring modified objective function
J = L - [0, 0, 1, 0]*mod;

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

%assign solver
options = struct;
options.ipopt.print_level = 0;
solver = nlpsol('solver','ipopt',nlp,options);

% Solve
sol = solver('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

% Extract Solution
w_opt_SS = full(sol.x);
%optimal inputs
uOptk = w_opt_SS(5);
%non modified prediction of the OF
costOFk = -w_opt_SS(3);
%modified prediction of the OF
costOFkMod = -full(sol.f);

end

