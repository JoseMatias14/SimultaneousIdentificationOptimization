function [lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas(par)

% variables
%1: choke outlet pressure % [Pa] -- 1e7
P_choke_out_lb = 1e-6;
P_choke_out_ub = 1e3;
    
%2: volumetric flow from the reservoir %[m3/s] -- 1e2
q_res_lb = 1e-6;
q_res_ub = 1e4;

%3: inlet volumetric flowrate % [m3/s] -- 1e2
q_comp_in_lb = q_res_lb;
q_comp_in_ub = q_res_ub;

%4: temperature % [K] -- 1e2
T_comp_out_lb = 1e-2;
T_comp_out_ub = 100;

%5: pressure % [Pa] -- 1e7
P_comp_out_lb = P_choke_out_lb;
P_comp_out_ub = 1e2*P_choke_out_ub;

%6: compressor efficiency
nu_lb = 1e-6;
nu_ub = 1.05;

%7: compressor head % [m] -- 1e3
H_lb = 1e-6;
H_ub = 1e4;

%8: compressor power %[W] - J/s
Pow_lb = 1e-12;
Pow_ub = 1e3;

%Pump
%9: inlet volumetric flowrate  % [m3/s] -- 1e2
q_pump_in_lb = q_res_lb;
q_pump_in_ub = q_res_ub;

%8: compressor power %[W] - J/s
Pop_lb = 1e-12;
Pop_ub = 1e3;

% inputs
%compressibility factor 
u_choke_lb = 0.5;
u_choke_ub = 1;

%gas lift rate [kg/s]
u_comp_lb = 0.75;
u_comp_ub = 1.05;


lbx = vertcat(P_choke_out_lb,q_res_lb,q_comp_in_lb,T_comp_out_lb,P_comp_out_lb,nu_lb,H_lb,Pow_lb,q_pump_in_lb,Pop_lb);
ubx = vertcat(P_choke_out_ub,q_res_ub,q_comp_in_ub,T_comp_out_ub,P_comp_out_ub,nu_ub,H_ub,Pow_ub,q_pump_in_ub,Pop_ub);
          
lbu = vertcat(u_choke_lb,u_comp_lb);
ubu = vertcat(u_choke_ub,u_comp_ub);
