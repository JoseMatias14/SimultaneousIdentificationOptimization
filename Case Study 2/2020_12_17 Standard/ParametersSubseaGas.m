function  par = ParametersSubseaGas
x = [.92,.05,.02,.005,.005,0,0,0,0,0,0];

% fixed
%choke valve constant
par.C_choke = 0.4*.4179402833086; %; [kg^0.5 m^0.5]
%separator diameter
par.D = 2; %[m]
%separator cross sectional area 
par.A = pi*(par.D/2)^2; %[m2]
%gas constant
par.R = 8.31446; % [J/(K mol)]
%gravity acceleration
par.g = 9.80665; % [m/s2]
%molar mass
% Molar masses [kg/kmol]
M_ =[...
    16.04, ... % C1
    30.07, ... % C2
    44.10, ... % C3
    58.12, ... % n-C4
    72.15, ... % n-C5
    86.18, ... % n-C6
    18.02, ... % H2O
    44.01, ... % CO2
    28.01, ... % N2
    58.12, ... % i-C4
    72.15, ... % i-C5
    ];

par.Mm = M_*x'; % [kg/kmol]
%liquid density - from the reservoir
par.rho_liq = 10; % [kg/m3]
%compressor head parameter
par.head_par = 4; % [m]

%pseudo critical mixture pressure
% Critical pressures [Pa]
Pc_ =[...
    4.60551724137931E6, ... % C1
    4.88137931034483E6, ... % C2
    4.25034482758621E6, ... % C3
    3.80000000000000E6, ... % n-C4
    3.37241379310345E6, ... % n-C5
    3.01379310344828E6, ... % n-C6
    22.0620689655172E6, ... % H2O
    7.38344827586207E6, ... % CO2
    3.39000000000000E6, ... % N2
    3.64896551724138E6, ... % i-C4
    3.39000000000000E6, ... % i-C5
    ];
par.Ppc = Pc_*x'; %[Pa]

%pseudo critical mixture temperature
% Critical temperatures [K]
Tc_ =[...
    190, ... % C1
    305, ... % C2
    370, ... % C3
    425, ... % n-C4
    469, ... % n-C5
    507, ... % n-C6
    647, ... % H2O
    304, ... % CO2
    126, ... % N2
    408, ... % i-C4
    460, ... % i-C5
    ];
par.Tpc = Tc_*x'; %[K]

%compressibility factor equation
par.a1 = 0.3265; % [-]
par.a2 = -1.0700; % [K]
par.a3 = -0.5339; % [K^3]
par.a4 = 0.01569; % [K^4]
par.a5 = -0.05165; % [K^5]
par.a6 = 0.5475; % [-]
par.a7 = -0.7361; % [K]
par.a8 = 0.1844; %[K^2]
par.a9 = 0.1056; % [-]
par.a10 = 0.6134; % [K^3]
par.a11 = 0.7210; % [-]

%Thermodynamic - cp calculation
coeffs =[...
    1.702,  9.081E-3,  -2.164E-6,       0; ... % C1
    1.131, 19.225E-3,  -5.561E-6,       0; ... % C2
    1.213, 28.785E-3,  -8.824E-6,       0; ... % C3
    1.935, 36.915E-3, -11.402E-6,       0; ... % n-C4
    2.464, 45.351E-3, -14.111E-6,       0; ... % n-C5
    3.025, 53.722E-3, -16.791E-6,       0; ... % n-C6
    3.470,  1.450E-3,          0, 0.121E5; ... % H2O
    5.457,  1.045E-3,          0,-1.157E5; ... % CO2
    3.280,  0.593E-3,          0, 0.040E5; ... % N2
    1.677, 37.853E-3, -11.945E-6,       0; ... % i-C4
    2.464, 45.351E-3, -14.111E-6,       0; ... % i-C5
    ];
c = x*coeffs;
par.c1 = c(1); %[-]
par.c2 = c(2); %[-/K]
par.c3 = c(3); %[-/K^2]
par.c4 = c(4); %[-/K^2]

%adjustable
%Compressor Head
par.c5 = -0.9937; % [s2/m6]
par.c6 = 2.256; % [s/m3]
par.c7 = 1.888; % [-]

%%%% 
% check these values
%%%%

% %Surge condition calculation
% par.c8 = 1; %[s2/m6]
% par.c9 = 1; %[s/m3]
% par.c10 = 1; %[-]
% 
% %Stonewall condition calculation
% par.c11 = 1; %[s2/m6]
% par.c12 = 1; %[s/m3]
% par.c13 = 1; %[-]

%separation efficiency
par.c14 = 5; %[kg/m3]
par.c15 = 7; %[(m^{15/2} s^3)/kg^{9/2}]

%compressor map constants
par.n11 = 0.582; % [s2/m6]
par.n12 = -2.398; % [s/m3]
par.n13 = 2.75; % [-]
par.n14 = -3.969; % [s/m3]
par.n15 = 4.303; % [-]
% par.n11 = 0.592; % [s2/m6]
% par.n12 = -2.398; % [s/m3]
% par.n13 = 2.75; % [-]
% par.n14 = -3.969; % [s/m3]
% par.n15 = 4.303; % [-]

%compressor map constants
par.n21 = -0.5; % []
par.n22 = 1.8; % []
par.n23 = -0.85; % []
% par.n21 = -0.55; % []
% par.n22 = 1.9; % []
% par.n23 = -0.85; % []

%compressor map constants
par.n31 = 0.43; % []
par.n32 = 0.7; % []
par.n33 = 5; % []
par.n34 = 1; % []
par.n35 = 2; % []
% par.n31 = 0.3; % []
% par.n32 = 0.9; % []
% par.n33 = 2; % []
% par.n34 = 1; % []
% par.n35 = 2; % []

%weights for the total modifier
% par.Q1 = diag([0,0,1e0,0,1e0,1e1,0,0]);
% par.Q2 = diag([0,0,1e0,0,1e0,1e1,0,0,0,0,1e0,0,1e0,1e1,0,0]);
par.Q1 = eye(8);
par.Q2 = eye(16);

%Optimization
par.qN_min = 1.163;
par.qN_max = 2.286;

par.Pcomp_min = 100*1e5; %100 bar

% changes in the plant model 
par.ph = 51;
par.ph_2_d1 = 131; %101
par.pd1 = 181; %151
par.pd1_2_d2 = 231; %201
par.pd2 = 250;

% Maintenance preiods
par.maintPeriods = 5;
par.stopPeriods = 25;








