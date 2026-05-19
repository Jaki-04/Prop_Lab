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
Pc = m_a * cp_c * (T_tot3 - T_tot2);
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
p1 = ptot_Out_diff / ((1 + ((g_GC - 1) * (M1 ^ 2)) / 2)) ^ (g_GC / (g_GC - 1));
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

%% Rotore

omega = 755.0118;   % da compressore
r_m = 0.4616;       % da compressore 
U_m = omega * r_m;

% Sezione d'ingresso
A2 = m_g / (rho2 * C2_z);
h2 = A2 / (2 * pi * r_m);
r_t2 = r_m + h2 / 2;
r_h2 = r_m - h2 / 2;

W2_z = C2_z;
W2_t = Ct_2 - U_m;
W2 = sqrt(W2_t^2 + W2_z^2);
a2 = sqrt(g_GC * R_GC * T2);

T2_t = T2 * (1 + (g_GC - 1) / 2 * M2 ^ 2);
p2_t = p2 * (1 + (g_GC - 1)/2 * M2^2)^(g_GC / (g_GC - 1));

% Sezione d'uscita
M_w3 = 0.8;
W3 = sqrt((a2^2 + (g_GC - 1) * (W2^2 / 2)) / ((1 / M_w3^2) + (g_GC - 1) / 2));

W3_z = C2_z;
W3_theta = sqrt(W3^2 - W3_z^2);

C3_t = U_m - W3_theta;
C3 = sqrt(C3_t^2 + W3_z^2);

R = 1 - (Ct_2 + C3_t) / (2 * U_m); % Grado di reazione

w_stadio = U_m * (Ct_2 - C3_t);
T_t3 = T2_t - (w_stadio / cp_GC);

% Temperatura statica alla sezione 3
T3 = T_t3 - (C3^2) / (2 * cp_GC);

% Numero di Mach in uscita
a3 = sqrt(g_GC * R_GC * T3);
M3 = C3 / a3;

% Pressioni alla sezione 3
p_t3 = p2_t * (T_t3 / T2_t)^(g_GC / ((g_GC - 1) * e_t));
p3 = p_t3 / (1 + (g_GC - 1)/2 * M3^2)^(g_GC / (g_GC - 1));

% Densità e Area di passaggio alla sezione 3
rho3 = p3 / (R_GC * T3);
A3 = m_g / (rho3 * W3_z);

% Geometria della pala in uscita
h3 = A3 / (2 * pi * r_m);
r_t3 = r_m + h3 / 2;
r_h3 = r_m - h3 / 2;

%% Sezione di scarico
T5_cool = T_t5cool - (0.5 * C2_z^2) / cp_GC;
M5 = C2_z / sqrt(g_GC * R_GC * T5_cool);
p5_tot = pi_t * ptot_Out_diff;
p5 = p5_tot / (1 + (M5^2) * (g_GC - 1) / 2);
rho5 = p5 / (R_GC * T5_cool);
A5 = m0 / (rho5 * C2_z);

h5 = A5 / (2 * pi * r_m);
r_t5 = r_m + h5 / 2;
r_h5 = r_m - h5 / 2;