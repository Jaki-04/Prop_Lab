clear
close
clc

%% Dimensionamento turbina

subCruise = Subcruise_parameters();
Air = Air_parameters('12000');

m_a = subCruise.m_a;
f = subCruise.f;
cp_GC = Air.cp_GC;
g_GC = Air.g_GC;
cp_c = Air.cp;
R_GC = Air.R_GC;

M2_sub = 0.5;
eta_m = 0.9;
eps_cool = 0.036;
m_g = m_a * (1 + f);
m_c = eps_cool * m_a;
m0 = m_g + m_c;
e_t = 0.9;

T_tot3 = 581.4252;
T_2 = 236.1485;
T_tot2 = T_2 * (1 + (g_GC - 1) / 2 * M2_sub ^ 2);
T_Out_cc = 1.3900e+03;
ptot_Out_diff = 4.3335e+05;

% Potenze compressore e turbina
Pc = m_a * cp_GC * (T_tot3 - T_tot2);
Pt = Pc / eta_m;

T_t5cool = ((1 - eps_cool) * cp_GC * T_Out_cc + eps_cool * cp_c * T_tot3 - Pt...
    / m0) / ((1 - eps_cool) * cp_GC + eps_cool * cp_c);
taut = T_t5cool / T_Out_cc;
pi_t = taut ^ (g_GC / ((g_GC - 1) * e_t));
SOT = 2 * T_Out_cc / (g_GC - 1);

%% Statore

M1 = M2_sub
T1 = T_Out_cc / (1 + (g_GC - 1) * (M1 ^ 2) / 2);
p1 = ptot_Out_diff / ((1 + (g_GC - 1) * (M1 ^ 2) / 2)) ^ (g_GC / (g_GC - 1));
C1 = M1 * sqrt(g_GC * R_GC * T1);