function  par = ParametersBlockReactor

%number of components
par.nc = 4;
%number of reaction sets available
par.nr = 9; % set 1 for only MAy

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
% ntot = 3;
% [a,b] = meshgrid(linspace(0.5,1.0,ntot),linspace(1.0,2.0,ntot));
% aflat = reshape(a,1,ntot^2);
% bflat = reshape(b,1,ntot^2);
% 
% for ii = 1:9
%     par.k{ii} = [aflat(ii);bflat(ii)];
% end

par.k = [0.75, 1.5,    0,  0;
    0.73, 1.45,    0,  0;
    0.75, 1.45,    0,  0;
    0,   0, 0.0254,  0.0166;
    0,   0,  0.02,  0.02;
    0,   0,  0.028,  0.02;
    0.1295,   0,    0,  0;
    0.135,   0,    0,  0;
    0.12,   0,    0,  0]';

%reaction constants
par.k_real = [0.75;1.5];

