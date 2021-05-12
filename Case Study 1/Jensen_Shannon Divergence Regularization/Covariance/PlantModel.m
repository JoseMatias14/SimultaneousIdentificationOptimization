function [diff,x_var,u_var,L,Fendxk,FendMeas,FendGradMeas] = PlantModel(dx0,u0,par)
%    Block reactor model - Implementation of the isothermal CSTR
%                Real reaction set: A + 2B -> C & B -> C (Model 1)
%
% Inputs:
%    dx0 = initial differential states
%    u0 = inputs
%    par = system parameters
%
% Outputs:
%    diff: system symbolic description (casadi)
%    x_var,u_var: symbolic states and inputs
%    L: symbolic quadrature value (Cc)
%    Fendxk: numerical differential states values
%    FendMeas: numerical mesurements
%    FendMeasNoise: numerical mesurements with noise
%    FendGradMeasu: numerical gradient of the mesurements
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br
% March 2018; Last revision: 13-Mar-2018
        
addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

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
R = vertcat(par.k_real(1)*x_var(1)*x_var(2),par.k_real(2)*x_var(2)^2);
 
%model
diff = vertcat(u_var*par.Cain/par.V,par.Fb*par.Cbin/par.V,0,0) - (u_var + par.Fb)/par.V*ones(4,1).*x_var + Upsilon*R;

%Objective function (optimization: maximize C)
L = -x_var(3); 

%end modeling

%% Casadi commands
%formalizing a function to be used in the *rootfinder*
residual = Function('residual',{x_var,u_var},{diff},{'x0','u'},{'xSS'});

%evaluating the roots of the system
%the order is: dummy label, resolution procedure, residual
Groot = rootfinder('Groot','kinsol',residual);% kinsol| nlpsol | newton

%now, we give the initial guess and initialize the rootfiner
sol = Groot('x0',dx0,'u',u0);

%extracting the results (from symbolic to numerical)
Fendxk = full(sol.xSS);
FendMeas = Fendxk + 0.02*randn(4,1); %Multiplying a random variable by a constant increases the variance by the square of the constant.

%assuming inputs as symbolic in order to obtain the gradients symbolically
%calling the integrator (symbolically)
Ffun = Groot('x0',dx0,'u',u_var);
xf = Ffun.xSS; %Css

%extracting SS gradient (symbolically)
gradMeas = Function('gradMeasu',{u_var},{jacobian(xf,u_var)});

%extracting the results (from symbolic to numerical)
FendGradMeas = full(gradMeas(u0)) + 0.02*randn(4,1);





