% Optimizing and obtaining the plant surface of the isothermal CSTR
% True reaction set: A + B -> C & 2B -> D

% Other m-files required: 
% InitialConditionBlockReactor.m; OptimizationBoundsBlockReactor.m;

% MAT-files required: none

clear 
close all
clc

%saving file
results = 'PlantSurface_BR';

%adding casadi
import casadi.*

%% Parameters
%number of components
par.nc = 4;
%number of reaction sets available
par.nr = 3;
%inflow of B
par.Fa = 5; %[L/min]
%inflow of B
par.Fb = 5; %[L/min]
%inlet concentration of A
par.Cain = 2; %[mol/L]
%inflow of B
par.Cbin = 1.5; %[mol/L]
%reactor volume
par.V = 500; %[L]
%reaction constants
par.k1 = 0.75; %[L/(mol min)^?]
%reaction constants
par.k2 = 1.5; %[L/(mol min)^?]

%% Model from PlantModel.m
%% Differential states
%symbolic declaration
%concentration array [Ca, Cb, Cc, Cd]
x_var = MX.sym('C',par.nc,1); % 1-4 [mol/L]

%% Input
%A inflow rate
u_var = MX.sym('Fa',1); %[L/min]

%% System equations
%stoichiometric coefficients
Upsilon = [-1,0;
           -1,-2;
            1,0;
            0,1];

%reaction array
R = vertcat(par.k1*x_var(1)*x_var(2),par.k2*x_var(2)^2);
 
%model
diff = vertcat(u_var*par.Cain/par.V,par.Fb*par.Cbin/par.V,0,0) - (u_var + par.Fb)/par.V*ones(4,1).*x_var + Upsilon*R;
 
%% Casadi commands
%formalizing a function to be used in the *rootfinder*
residual = Function('residual',{x_var,u_var},{diff},{'x0','u'},{'xSS'});

%evaluating the zeros of the system
%rootfinder input order is: dummy label, resolution procedure, residual
Groot = rootfinder('Groot','kinsol',residual);% kinsol| nlpsol | newton

%% preparing a grid to see the dependence of Cc_SS to Fa
ntot = 200; %different Fa
nu1grid = linspace(0.05,10,ntot);

%creating map to compute the rootfinder automatically instead of using a loop
Grootgridflat = Groot.map(ntot);

% map: dummy initial guess for the states & fixed grid for u
out = Grootgridflat('x0',ones(par.nc,ntot),'u',nu1grid);

% extract the value of the OF
grid = full(out.xSS(3,:));

%Plotting
figure(1)

plot(nu1grid,grid)
hold on 

title('Dependence of C_{SS} to Fa')
xlabel('Inlet of A [L/min]')
ylabel('Concentration of C [mol/L]')


%% Optimization 
%states and inputs bounds
[lbx,lbu,ubx,ubu] = OptimizationBoundsBlockReactor;
%initial condition
[dx0,u0] = InitialConditionBlockReactor;

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
lbw = [lbw;lbx;lbu];
ubw = [ubw;ubx;ubu];
wk = [wk;dx0;u0];
% Note: Equality constraint is defined by setting ub = lb with the same
% value

%Add the system model as constraints
g = {g{:},diff};
lbg = [lbg;zeros(par.nc,1)]; %SS -> dif and alg == 0
ubg = [ubg;zeros(par.nc,1)];

L = -x_var(3);%maximize Cc

% formalize it into an NLP problem
nlp = struct('x',vertcat(w{:}),'f',L,'g',vertcat(g{:}));

% Assign solver
solver = nlpsol('solver','ipopt',nlp);

% Solve
sol = solver('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

% Extract Solution
w_opt_SS = full(sol.x);
uOptk = w_opt_SS(5);
cost_OF = -full(sol.f);

plot(uOptk,cost_OF,'r*')

%Saving results
save(results,'cost_OF','nu1grid','grid','uOptk');


