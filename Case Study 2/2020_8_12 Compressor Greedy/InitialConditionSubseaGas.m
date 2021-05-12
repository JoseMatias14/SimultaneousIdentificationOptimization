function [x0,u0,f0] = InitialConditionSubseaGas(par)

%% 
% parameters
% x = [.92,.05,.02,.005,.005,0,0,0,0,0,0];
% par = ParametersSubseaGas;

%inputs
%choke valve openings 
u_choke = 0.8; % u_1 = 0.565 [%]
%compressor speed (normalized)
u_comp = 0.9; % u_2 = 0.85[%]

u0 = vertcat(u_choke,u_comp);

%feed
%mass flow gas from the reservoir
m_res_gas = 0.9; %0.9|0.8997 [kg/s] -- 1e2
%mass flow liquid from the reservoir
m_res_liq = 0.3; %0|1.4616[kg/s] -- 1e2
%reservoir pressure
P_res = 1.1; % 1.0 [Pa]  -- 1e7
%reservoir temperature
T_res = 3.5; % 3.5 [K]  -- 1e2 

f0 = vertcat(m_res_gas,m_res_liq,P_res,T_res);

% States
% %choke outlet pressure
% --> P_choke_out [Pa] -- 1e7
% %volumetric flow from the reservoir
% --> q_res [m3/s]
% %Compressor
% %inlet volumetric flowrate
% -->  q_comp_in [m3/s]
% %temperature
%  --> T_comp_out [K] -- 1e2
% %pressure
%  --> P_comp_out [Pa] -- 1e7
% %compressor efficiency
%  --> nu [-]
% %compressor head
%  --> H  [m] -- 1e1
% %compressor power
%  --> Pow [W] -- 1e4
% %Pump
% %inlet volumetric flowrate
% q_pump_in [m3/s] -- 1e-1
% %compressor power
% Pop [W] -- 1e3
xGuess = [0.987913152093562;1.49990184193191;1.46159550354337;3.49960000000000;0.988360495898890;0.798712028891691;0.815176284684305;0.883677989303943;0.000383063383885381;1.71361031829183e-05];
[x0,~] = PlantModel(xGuess,u0,f0,par,0);

%clear IPOPT output
clc
end

