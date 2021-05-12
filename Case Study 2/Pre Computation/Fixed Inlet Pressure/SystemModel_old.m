function [yk,grad_yk,alg,x_var,u_var,f_var] = SystemModel(xk_1,uk,fk,par,modelFlag,thetaHatk)
% Build system model for optimization (IPOPT)
% Subsea gas compression model - The purpose of the gas compression station is to
%                               boost the pressure of the stream so that it is sufficiently
%                               high to overcome the pressure drop in the transportation
%                               pipeline and arrive at the receiving facility topside with
%                               the desired outlet pressure.
%                               based on: Health Estimation and Optimal
%                               Operation of a Subsea Gas Compression
%                               System Under Uncertainty - A. Verheyleweghen
%
% Inputs:
%|      par = system parameters
%
% Outputs:
%       alg = steady-state system model (CASADI)
%       x_var = system states (CASADI)
%       u_var = system inputs (CASADI)
%       f_var = system feed (CASADI)

%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Jose Otavio Matias
% email: jose.otavio@usp.br
% August 2019; Last revision: 08-Aug-2019

addpath ('\\home.ansatt.ntnu.no\joseoa\Documents\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

    %% Parameters
    %x = [.92,.05,.02,.005,.005,0,0,0,0,0,0];
    % fixed
    %choke valve constant
    C_choke = par.C_choke; %0.4*.4179402833086; [kg^0.5 m^0.5]
    %separator cross sectional area 
    A = pi*(par.D/2)^2; %[m2]
    %gas constant
    R = par.R; %8.31446 [J/(K mol)]
    %gravity acceleration
    g = par.g; %9.80665 [m/s2]
    %molar mass
    Mm = par.Mm; % [g/mol]
    %liquid density - from the reservoir
    rho_liq = par.rho_liq ; %10 [kg/m3]
    %compressor head parameter
    head_par = par.head_par; %4 [?????]

%     %pseudo critical mixture pressure
%     Ppc = par.Ppc; %[Pa]
%     %pseudo critical mixture temperature
%     Tpc = par.Tpc; %[K]
%     %compressibility factor equation
%     a1 = par.a1; %0.3265 [-]
%     a2 = par.a2; %-1.0700 [K]
%     a3 = par.a3; %-0.5339[K^3]
%     a4 = par.a4; %0.01569 [K^4]
%     a5 = par.a5; %-0.05165 [K^5]
%     a6 = par.a6; %0.5475 [-]
%     a7 = par.a7; %-0.7361 [K]
%     a8 = par.a8; %0.1844[K^2]
%     a9 = par.a9; %0.1056 [-]
%     a10 = par.a10; %0.6134 [K^3]
%     a11 = par.a11; %0.7210 [-]

    %Thermodynamic - cp calculation
    c1 = par.c1; %[-]
    c2 = par.c2; %[1/K]
    c3 = par.c3; %[1/K^2]
    c4 = par.c4; %[K^2]

    %adjustable
    %Compressor Head
    c5 = SX.sym('c5'); %-0.9937 [s2/m6]
    c6 = SX.sym('c6'); %2.256 [s/m3]
    c7 = SX.sym('c7'); %1.888 [-]
%     %Surge condition calculation
%     c8 = SX.sym('c8'); %[s2/m6]
%     c9 = SX.sym('c9'); %[s/m3]
%     c10 = SX.sym('c10'); %[-]
%     %Stonewall condition calculation
%     c11 = SX.sym('c11'); %[s2/m6]
%     c12 = SX.sym('c12'); %[s/m3]
%     c13 = SX.sym('c13'); %[-]
    %separation efficiency
    c14 = SX.sym('c14'); %5 [kg2/m6]
    c15 = SX.sym('c15'); %7 [(m^{13/2} s)/kg^{5/2}]
    %compressor map constants
    n11 = SX.sym('n11'); 
    n12 = SX.sym('n12'); 
    n13 = SX.sym('n13'); 
    n14 = SX.sym('n14'); 
    n15 = SX.sym('n15'); 
    
    n21 = SX.sym('n21'); 
    n22 = SX.sym('n22'); 
    n23 = SX.sym('n23'); 
    
    n31 = SX.sym('n31'); 
    n32 = SX.sym('n32'); 
    n33 = SX.sym('n33'); 
    n34 = SX.sym('n34'); 
    n35 = SX.sym('n35'); 

    %% Inputs
    %Manipulated variables
    %choke valve openings 
    u_choke = SX.sym('u_choke'); % u_1 = 0.565 [%]
    %compressor speed (normalized)
    u_comp = SX.sym('u_comp'); % u_2 = 0.85[%]

    % Feed
    %mass flow gas from the reservoir
    m_res_gas = SX.sym('m_res_gas'); %0|1.4616 [kg/s] -- 1e2
    %mass flow liquid from the reservoir
    m_res_liq = SX.sym('m_res_liq'); %0.9|0.8997 [kg/s] -- 1e2
    %reservoir pressure
    P_res = SX.sym('P_res'); % 1.0 [Pa]  -- 1e7
    %reservoir temperature
    T_res = SX.sym('T_res'); % 3.5 [K]  -- 1e2 

    %% Variables
    %Reservoir
    %mass flow from the reservoir
    %volumetric flow from the reservoir
    q_res = SX.sym('q_res'); %[m3/s]

    %Choke
    %choke outlet pressure
    P_choke_out = SX.sym('P_choke_out'); % [Pa] -- 1e7

    %Compressor
    %inlet volumetric flowrate
    q_comp_in = SX.sym('q_comp_in'); % [m3/s]
    %temperature
    T_comp_out = SX.sym('T_comp_out'); % [K] -- 1e2
    %pressure
    P_comp_out = SX.sym('P_comp_out'); % [Pa] -- 1e7
    %compressor efficiency
    nu = SX.sym('nu'); %[-]
    %compressor head
    H = SX.sym('H'); % [m] -- 1e1
    %compressor power
    Pow = SX.sym('Pow'); %[W] -- 1e4

    %Pump
    %inlet volumetric flowrate
    q_pump_in = SX.sym('q_pump_in'); % [m3/s] -- 1e-1
    %compressor power
    Pop = SX.sym('Pop'); %[W] -- 1e3
    
    %% Modeling
    %supporting relations
    % ===================================
    %     Choke
    % ===================================
    %Density %[kg/m3] -- 1e2
    rho_gas_choke_out = 1e-2*(P_choke_out*1e7)*(Mm*1e-3)/(8.3144*(T_res*1e2));

    % ===================================
    %     Separator
    % ===================================
    %separator efficiency K  %[kg^0.5/(s m^0.5)]
    K = (q_res)/A*sqrt((rho_gas_choke_out*1e2)/((rho_liq*1e2) - (rho_gas_choke_out*1e2)));
    %separator efficiency - related to mass %[-]
    alpha = (1 - (rho_gas_choke_out/c14)^2 - c15*rho_gas_choke_out^2*K^3);
    %heat capaticy of the outlet stream %[J/(K mol)]
    Cp_sep = (c1 + c2*(T_res*1e2) + c3*(T_res*1e2)^2 + c4*(T_res*1e2)^(-2))*R;
    % ===================================
    %     Compressor
    % ===================================
    %gas volume fraction %[-]
    GVF = (alpha*(m_res_gas*1e2)/(rho_gas_choke_out*1e2))/(q_comp_in);%qgas/(qliq + qgas)
    %compressor average density at the inlet %[kg/m3] -- 1e2
    rho_comp_avg = 1e-2*(GVF*(rho_gas_choke_out*1e2) + (1 - GVF)*(rho_liq*1e2));
    %heat capaticy of the outlet stream %[J/(K mol)] 
    Cp_comp = (c1 + c2*(T_comp_out*1e2) + c3*(T_comp_out*1e2)^2 + c4*(T_comp_out*1e2)^(-2))*R;
    %k is defined in terms of adiabatic ratios times the efficiency %[-]
    gamma = 0.5*(Cp_sep/(Cp_sep - R) + Cp_comp/(Cp_comp - R));
    k = nu*gamma/(1 - gamma);
    % Flow corresponding to N=1 (fan laws) % [m3/s]
    qN = (q_comp_in)/u_comp;
    %wood correction factor %[-]
    fwood  = 1/((rho_comp_avg*1e2)/(rho_gas_choke_out*1e2)*sqrt(GVF*(rho_comp_avg*1e2)/(rho_gas_choke_out*1e2)));

    %Steady-state model
    %choke valve characteristics
    f1 = ((m_res_gas*1e2) + (m_res_liq*1e2)) - u_choke*C_choke*sqrt((P_res*1e7) - (P_choke_out*1e7));
    %volumetric flowrates
    f2 = (q_res) - ((m_res_gas*1e2)/(rho_gas_choke_out*1e2) + (m_res_liq*1e2)/(rho_liq*1e2));
    %compressor inlet volumetric flowrate
    f3 = (q_comp_in) - (alpha*(m_res_gas*1e2)/(rho_gas_choke_out*1e2) + (1 - alpha)*(m_res_liq*1e2)/(rho_liq*1e2));
    %polytropic relation
    f4 = P_comp_out - P_choke_out*(T_comp_out/T_res)^k;
    %compressor efficiency - calculated by the compressor map
    if modelFlag == 1
        f5 = nu - ((n11*qN^2 + n12*qN + n13)/(qN^2 + n14*qN + n15));
    elseif modelFlag == 2
        f5 = nu - (n21*qN^2 + n22*qN + n23);
    else %modelFlag == 3
        f5 = nu - (n31 + n32*qN.^2./(n33 + qN.^2.*(n34 + qN.^2/n35)));
    end
    
    %compressor head (source ?)
    f6 = (H*1e1) - head_par*(c5*qN^2 + c6*qN + c7)*u_comp^2*fwood;
    %compressor power
    f7 = (Pow*1e4) - (H*1e1)*(q_comp_in)*(rho_gas_choke_out*1e2)*g/nu;
    %temperature outlet (source ?)
    f8 = (H*1e1) - k*R/(g*(Mm*1e-3))*((T_comp_out*1e2) - (T_res*1e2));
    %pump inlet volumetric flowrate
    f9 = (q_pump_in*1e-1) - ((1 - alpha)*(m_res_gas*1e2)/(rho_gas_choke_out*1e2) + alpha*(m_res_liq*1e2)/(rho_liq*1e2));
    %pump power
    f10 = (Pop*1e3) - ((q_pump_in*1e-1)*((P_comp_out*1e7) - (P_choke_out*1e7))); %[W] -- 1e7

    % Form the equation system
    alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10);

    % give parameter values
    alg = substitute(alg,c5,par.c5);
    alg = substitute(alg,c6,par.c6);
    alg = substitute(alg,c7,par.c7);
%     alg = substitute(alg,c8,par.c8);
%     alg = substitute(alg,c9,par.c9);
%     alg = substitute(alg,c10,par.c10);
%     alg = substitute(alg,c11,par.c11);
%     alg = substitute(alg,c12,par.c12);
%     alg = substitute(alg,c13,par.c13);
    alg = substitute(alg,c14,par.c14);
    alg = substitute(alg,c15,par.c15);
    if modelFlag == 1
        alg = substitute(alg,n11,thetaHatk(1));
        alg = substitute(alg,n12,thetaHatk(2));
        alg = substitute(alg,n13,thetaHatk(3));
        alg = substitute(alg,n14,thetaHatk(4));
        alg = substitute(alg,n15,thetaHatk(5));
    elseif modelFlag == 2
        alg = substitute(alg,n21,thetaHatk(1));
        alg = substitute(alg,n22,thetaHatk(2));
        alg = substitute(alg,n23,thetaHatk(3));
    else %modelFlag == 3
        alg = substitute(alg,n31,thetaHatk(1));
        alg = substitute(alg,n32,thetaHatk(2));
        alg = substitute(alg,n33,thetaHatk(3));
        alg = substitute(alg,n34,thetaHatk(4));
        alg = substitute(alg,n35,thetaHatk(5));
    end

    % concatenate the differential and algebraic states
    x_var = vertcat(P_choke_out,q_res,q_comp_in,T_comp_out,P_comp_out,nu,H,Pow,q_pump_in,Pop); %[-]
    u_var = vertcat(u_choke,u_comp);
    f_var = vertcat(m_res_gas,m_res_liq,P_res,T_res);

    %end modeling

    %% Casadi commands
    % Building NLP
    % decision variables
    w = {};

    %first input corresponds to the empty cell
    % second to the system variables
    % third to the system inputs
    % fourth to the system feed variables
    w = {w{:},x_var,u_var,f_var};

    % constraints
    g = {};
    %Add the system model as constraints
    g = {g{:},alg};

    % objective function
    J = 0;

    % formalize it into an NLP problem
    nlp = struct('x',vertcat(w{:}),'f',J,'g',vertcat(g{:}));

    % Assign solver
    options = struct;
    options.ipopt.print_level = 5;
%    options.ipopt.return_status = true;
    F = nlpsol('solver','ipopt',nlp,options);

    % ===================================
    %     SS simulation
    % ===================================
    % Solving
    %guess for the varibles
    %[xGuess,~,~] = InitialGuessSubseaGasCalculation(par);
    %setting lower and upper bounds
    [lbx,ubx,~,~] = OptimizationBoundsSubseaGas(par);

    %specifying values
    wk = [];
    lbw = [];
    ubw = [];

    lbw = [lbw;lbx;uk;fk];%one is a dummy variable - represent time
    ubw = [ubw;ubx;uk;fk];
    wk = [wk;xk_1;uk;fk];

    lbg = [];
    ubg = [];
    lbg = [lbg;zeros(10,1)]; %SS - dif and alg == 0
    ubg = [ubg;zeros(10,1)];

    % Solve
    sol = F('x0',wk,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);

    clc
    % Extract Solution
    xSol = full(sol.x);
    %extracting the results (from symbolic to numerical)
    xk = xSol(1:10);
    %model predictions
    yk = par.H*xk;
    
    
    % ===================================
    %     Gradient
    % ===================================
    %calculating gradients
    sens = Function('grad_uk',{x_var,u_var,f_var},{-jacobian(alg,x_var)\jacobian(alg,u_var)});

    % ===================================
    %     Gradient
    % ===================================
    grad_uk = full(sens(xk,uk,fk));
    
    %model gradients
    grad_yk = par.H*grad_uk;
    
end
