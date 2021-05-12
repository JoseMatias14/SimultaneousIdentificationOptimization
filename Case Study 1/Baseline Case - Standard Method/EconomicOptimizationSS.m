function [uOptk,costOFk,costOFkMod] = EconomicOptimizationSS(dx0,u0,par,diff,x_var,u_var,L,mod)

% Optmize the steady state model based on an economic criteria

% Inputs:
%    dx0,u0,par - initial guess for the states and inputs
%    par - system parameters
%    diff - symbolic model (differential equations)
%    x_var,u_var - symbolic states and inputs
%    L - cost function
%    mod - modifier
%
% Outputs:
%    uOptk - optimal inputs
%    costOFk - predicted value of the OF using the non-modified model
%    costOFkMod - objective function (cost) - modified model
%
% Subfunctions: OptimizationBoundsBlockReactor

% Author: Jose Matias
% email: jose.otavio@usp.br
% March 2018; Last revision: 13-Mar-2018

addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
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

w = {w{:},x_var,u_var};%first input corresponds to the empty cell
lbw = [lbw;-inf;-inf;-inf;-inf;lbu];
ubw = [ubw;inf;inf;inf;inf;ubu];
wk = [wk;dx0;u0];
% Note: Equality constraint is defined by setting ub = lb with the same
% value

%Add the system model as constraints
g = {g{:},diff};
lbg = [lbg;zeros(par.nc,1)]; %SS -> dif and alg == 0
ubg = [ubg;zeros(par.nc,1)];
% Note: Equality constraint is defined by setting ub = lb with the same
% value

%bounds on the modified states
g = {g{:},x_var + mod};
lbg = [lbg;lbx]; 
ubg = [ubg;ubx];

% formalize it into an NLP problem
J = L - [0, 0, 1, 0]*mod;
nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

%assign solver
options = struct;
options.ipopt.print_level = 0;
solver = nlpsol('solver','ipopt',nlp,options);

% Solve
sol = solver('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

% % nlp_hess_l:(x[118],p[],lam_f,lam_g[114])->(hess_gamma_x_x[118x118,1572nz]) MXFunction
% % ===================================
% %     Hessian
% % ===================================
% hf = solver.get_function('nlp_hess_l');
% %nlp_f,nlp_g,nlp_grad,nlp_grad_f,nlp_hess_l,nlp_jac_g.
% HessX = full(hf(full(sol.x),[],[],full(sol.lam_g)));
% 
% Juu = HessX(5,5);
% 
% gf = solver.get_function('nlp_grad_f');
% gradFull = full(gf(full(sol.x),[]));
% 
% %hf = Function('gradMeasu',{x_var,z_var,p_var,grad1,grad2},{J.hessian(GOR)});
% 
% lam_g = sol.lam_g;
% Lag = J + lam_g'*vertcat(g{:});
% % 
% gl = Function('gl',{x_var,u_var},{jacobian(Lag,u_var)});
% gLag = full(gl(full(sol.x(1:4)),full(sol.x(5))));
% % 
% hl = Function('gl',{x_var,u_var},{hessian(Lag,u_var)});
% hLag = full(hl(full(sol.x(1:4)),full(sol.x(5))));
% 
% % hf = Function('gradMeasu',{x_var,z_var,p_var,grad1,grad2},{hessian(dot(jacobian(Lag,p_res),jacobian(Lag,p_res)),p_res)});
% % Hess = full(hf(xEst(1:8),xEst(9:38),[uk;thetaHat],xEst(43:80),xEst(81:end)));


% Extract Solution
w_opt_SS = full(sol.x);
%optimal inputs
uOptk = w_opt_SS(5);
%non modified prediction of the OF
costOFk = -w_opt_SS(3);
%modified prediction of the OF
costOFkMod = -full(sol.f);

end

