clear
close
clc

%% Dimensionamento turbina

supCruise = Supercruise_parameters();
Air = Air_parameters('12000');

m_a = supCruise.m_a;
cp_GC = Air.cp_GC;
g_GC = Air.g_GC;
cp_c = Air.cp;

M2_sub = 0.5;
eta_m = 0.9;
eps_cool = 0.036;
m_0 = 

T_tot3 = 581.4252;
T_2 = 236.1485;
T_tot2 = T_2 * (1 + (g_GC - 1) / 2 * M2_sub ^ 2);
T_Out_cc = 1.3900e+03;

% Potenze compressore e turbina
Pc = m_a * cp_GC * (T_tot3 - T_tot2);
Pt = Pc / eta_m;
T_t5cool = ((1 - eps_cool) * cp_GC * T_Out_cc + eps_cool * cp_c * T_tot3 - Pt...
    / m_a) / ((1 - eps_cool) * cp_GC + eps_cool * cp_c);