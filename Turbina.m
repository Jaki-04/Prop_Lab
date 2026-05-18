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
eps_cool = 0.05;
m_g = m_a * (1 + f);
m_c = eps_cool * m_a;
m0 = m_g + m_c;
e_t = 0.9;

T_tot3 = 706.7317;
T_2 = 247.9559;
T_tot2 = T_2 * (1 + (g_GC - 1) / 2 * M2_sub ^ 2);
T_Out_cc = 1.7653e+03;
T_Out_cc_tot = T_Out_cc * (1 + (g_GC - 1) / 2 * M2_sub ^ 2);
ptot_Out_diff = 7.7935e+05;

% Potenze compressore e turbina
Pc = m_a * cp_GC * (T_tot3 - T_tot2);
Pt = Pc / eta_m;

T_t5cool = ((1 - eps_cool) * cp_GC * T_Out_cc_tot + eps_cool * cp_c * T_tot3 - Pt...
    / m0) / ((1 - eps_cool) * cp_GC + eps_cool * cp_c);
taut = T_t5cool / T_Out_cc_tot;
pi_t = taut ^ (g_GC / ((g_GC - 1) * e_t));
SOT = 2 * T_Out_cc_tot / (g_GC - 1);

%% Statore

% Ingresso statore
M1 = M2_sub;
T1 = T_Out_cc_tot / (1 + ((g_GC - 1) * (M1 ^ 2)) / 2);
p1 = ptot_Out_diff / (((1 + (g_GC - 1) * (M1 ^ 2)) / 2)) ^ (g_GC / (g_GC - 1));
C1 = M1 * sqrt(g_GC * R_GC * T1);
rho1 = p1 / (R_GC * T1);
A1 = m_g / (rho1 * C1);

M2 = 1.0;          
C2_z = C1;         

% Termodinamica in uscita
T2 = T_Out_cc_tot * (2 / (g_GC + 1));
p2 = ptot_Out_diff * (2 / (g_GC + 1)) ^ (g_GC / (g_GC - 1));
rho2 = p2 / (R_GC * T2);

C2 = sqrt(g_GC * R_GC * T2);
Ct_2 = sqrt(C2^2 - C2_z^2); % Velocità tangenziale

alpha2 = atan(Ct_2 / C2_z);
alpha2_deg = rad2deg(alpha2);

A_g = m_g / (rho2 * C2);    
