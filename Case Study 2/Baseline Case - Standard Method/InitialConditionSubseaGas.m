function [x0,u0,f0] = InitialConditionSubseaGas(par)
%   Initial condition for simulating the Subsea Compression Station model
%
% Inputs:
%    par: system parameters
% Outputs:
%    u0: initial input value
%    f0: initial feed value
%    x0: steady-state system states related to u0

% Subfunctions: PlantModel.m

%inputs
%choke valve openings 
u_choke = 0.8; %[-]
%compressor speed (normalized)
u_comp = 0.9; %[-]

u0 = vertcat(u_choke,u_comp);

%feed
%mass flow gas from the reservoir
m_res_gas = 0.9; %[1e2 kg/s]
%mass flow liquid from the reservoir
m_res_liq = 0.3; %[1e2 kg/s]
%reservoir pressure
P_res = 1.1; % [1e7 Pa]
%reservoir temperature
T_res = 3.5; % [1e2 K] 

f0 = vertcat(m_res_gas,m_res_liq,P_res,T_res);

% States
% %choke outlet pressure
% --> P_choke_out [1e7 Pa]
% %volumetric flow from the reservoir
% --> q_res [m3/s]
% %Compressor
% %inlet volumetric flowrate
% -->  q_comp_in [m3/s]
% %temperature
%  --> T_comp_out [1e2 K]
% %pressure
%  --> P_comp_out [1e7 Pa]
% %compressor efficiency
%  --> nu [-]
% %compressor head
%  --> H  [1e1 m]
% %compressor power
%  --> Pow [1e4 W]
% %Pump
% %inlet volumetric flowrate
% q_pump_in [1e-1 m3/s]
% %compressor power
% Pop [1e3 W]
xGuess = [0.987913152093562;1.49990184193191;1.46159550354337;3.49960000000000;0.988360495898890;0.798712028891691;0.815176284684305;0.883677989303943;0.000383063383885381;1.71361031829183e-05];
[x0,~] = PlantModel(xGuess,u0,f0,par,0);

%clear IPOPT output
clc
end

