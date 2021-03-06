function [lbx,lbu,ubx,ubu] = OptimizationBoundsBlockReactor
%   Bounds on the variables of the CSTR model (ReactorModel.m)
%
% Inputs:
%
% Outputs:
%    lbx, ubx: state bounds
%    lbu, ubu: input bounds

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

%states: concentration [Ca, Cb, Cc, Cd]
Ca_ub = 1e1; %[mol/L]
Cb_ub = 1e1; %[mol/L]
Cc_ub = 1e1; %[mol/L]
Cd_ub = 1e1; %[mol/L]

%A Inflow
Fa_ub = 10; %[L/min]

    ubx = vertcat(Ca_ub,Cb_ub,Cc_ub,Cd_ub);
    ubu = Fa_ub;

%states: concentration [Ca, Cb, Cc, Cd]
Ca_lb = 1e-8; %[mol/L]
Cb_lb = 1e-8; %[mol/L]
Cc_lb = 1e-8; %[mol/L]
Cd_lb = 1e-8; %[mol/L]

%A Inflow
Fa_lb = 0.1; %[L/min]
    
    lbx = vertcat(Ca_lb,Cb_lb,Cc_lb,Cd_lb);
    lbu = Fa_lb;


