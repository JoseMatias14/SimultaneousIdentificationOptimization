function [dx0,u0] = InitialConditionBlockReactor
%   Initial condition for simulating CSTR model
%
% Inputs:
%
% Outputs:
%    u0: initial input value
%    dx0: steady-state system states related to u0

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

%concentration array [Ca, Cb, Cc, Cd]
Ca0 = 0.892218722807789;    %[mol/L]
Cb0 = 0.0197760750473837;   %[mol/L]
Cc0 = 0.441114610525545;    %[mol/L]
Cd0 = 0.0195546572135357;   %[mol/L] 

%A Inflow
Fa0 = 10;                   %[L/min] 

dx0 = vertcat(Ca0,Cb0,Cc0,Cd0);
u0 = Fa0;