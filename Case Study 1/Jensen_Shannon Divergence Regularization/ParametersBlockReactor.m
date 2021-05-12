function  par = ParametersBlockReactor

%number of components
par.nc = 4;
%number of reaction sets available
par.nr = 9;

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
%reaction constants (for models)
par.k = [0.75, 1.2,    0,  0;
    0.72, 1.5,    0,  0;
    %0.75, 1.45,    0,  0;
    0.75,  1.5,    0,  0;
    0,   0, 0.0254,  0.0166;
    0,   0,  0.02,  0.02;
    0,   0,  0.028,  0.02;
    0.1295,   0,    0,  0;
    0.135,   0,    0,  0;
    0.12,   0,    0,  0]';

%reaction constants (for plant)
par.k_real = [0.75;1.5];

%number of points for computing the expected value of the objective function (Eq. 12) 
par.nint = 10;