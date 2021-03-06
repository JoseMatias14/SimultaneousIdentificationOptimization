function [diff,x_var,u_var,L,Fendxk,FendMeas,FendGradMeas] = ReactorModel(dx0,u0,par,flag)
%    Block reactor model - Implementation of the isothermal CSTR
%
% Inputs:
%    dx0 = initial differential states (guess)
%    u0 = inputs (fixed for simulation)
%    par = system parameters
%    flag = indicates which model structure to use
%
% Outputs:
%    diff: system symbolic description (casadi)
%    x_var,u_var: symbolic states and inputs
%    L: symbolic quadrature value (-Cc)
%    Fendxk: numerical differential states values
%    FendMeas: numerical mesurements
%    FendGradMeasu: numerical gradient
%
% Other m-files required: none
% MAT-files required: none
        
import casadi.*

% Multimodel tuning
% m1, m2, m3: A + B -k1-> C & 2B -k2-> D & pure D
% m1: k1 = 0.75,  k2 = 1.2 
% m2: k1 = 0.72,  k2 = 1.5 
% m3: k1 = 0.75,  k2 = 1.5  (TRUE MODEL)

% m4, m5, m6: A + B -k3-> 2C & C + 2B -k4-> D & pure D
% m4: k3 = 0.0254, k4 = 0.0166 
% m5: k3 = 0.02,   k4 = 0.02 
% m6: k3 = 0.028,  k4 = 0.02 

% m7, m8, m9: A + B -k1-> C & D as impurity
% m1: k1 = 0.1295
% m2: k1 = 0.135
% m3: k1 = 0.12 

%% Differential states
%symbolic declaration
%concentration array [Ca, Cb, Cc, Cd]
x_var = MX.sym('C',par.nc,1); % 1-4 [mol/L]

%% Input
%A inflow rate
u_var = MX.sym('Fa',1); %[L/min]
%Model determination (boolean)
xB = MX.sym('xB',par.nr); % 1-3 [-]

%% System equations
%%stoichiometric coefficients (adapt for the multimodel structure)
r1 = xB(1) + xB(2) + xB(3);
r2 = xB(4) + xB(5) + xB(6);
r3 = xB(7) + xB(8) + xB(9);

Upsilon = [-(r1 + r3),    0,  -r2,     0;
           -(r1 + r3),-2*r1,  -r2, -2*r2;
            (r1 + r3),    0, 2*r2,   -r2;
                    0,   r1,    0,    r2];
        
%reaction array
k1 = 0;
k2 = 0;
k3 = 0;
k4 = 0;
for ii = 1:par.nr
    k1 = k1 + xB(ii)*par.k(1,ii);
    k2 = k2 + xB(ii)*par.k(2,ii);
    k3 = k3 + xB(ii)*par.k(3,ii);
    k4 = k4 + xB(ii)*par.k(4,ii);
end

R = vertcat(k1*x_var(1)*x_var(2),k2*x_var(2)^2,k3*x_var(1)*x_var(2)^2,k4*x_var(3)*x_var(2)^2);

%model
diff = vertcat(u_var*par.Cain/par.V,par.Fb*par.Cbin/par.V,0,(r3*par.Fb*dx0(4)/par.V)) - (u_var + par.Fb)/par.V.*x_var + Upsilon*R;
       
%Objective function (optimization: maximize C)
L = -x_var(3); 
%end modeling

% give parameter values -> Determining the model structure
diff = substitute(diff,xB,flag);

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
FendMeas = Fendxk;

%assuming inputs as symbolic in order to obtain the gradients symbolically
%calling the integrator (symbolically)
Ffun = Groot('x0',dx0,'u',u_var);
xf = Ffun.xSS; %Css

%extracting SS gradient (symbolically)
gradMeas = Function('gradMeasu',{u_var},{jacobian(xf,u_var)});

%extracting the results (from symbolic to numerical)
FendGradMeas = full(gradMeas(u0));





