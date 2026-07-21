function [p_tOut, T_tOut, M_Out, A_Out_Comp, L_turbina] = Turbina()

subCruise = Subcruise_parameters();
Air = Air_parameters('12000');

m_a = subCruise.m_a;
f = subCruise.f;
cp_GC = Air.cp_GC;
g_GC = Air.g_GC;
cp_a = Air.cp;
R_GC = Air.R_GC;

M4 = 0.5;
eta_m = 0.92;
eps_cool = 0.0245;
m_0 = m_a * (1 + f);
e_t = 0.9;

T_t4=1420;
P_t4=4.17e+05;

%% Statore

% Ingresso statore
T4 = T_t4 / (1 + ((g_GC - 1) * (M4 ^ 2)) / 2);
P4 = P_t4 / ((1 + ((g_GC - 1) * (M4 ^ 2)) / 2)) ^ (g_GC / (g_GC - 1));
C1 = M4 * sqrt(g_GC * R_GC * T4);
rho4 = P4 / (R_GC * T4);
A1 = m_0 / (rho4 * C1);

Mstar = 1.0;          
Cz_2 = C1;         

% Termodinamica in uscita
T2 = T_t4* (2 / (g_GC + 1));
P2 = P_t4 * (2 / (g_GC + 1)) ^ (g_GC / (g_GC - 1));
rho2 = P2 / (R_GC * T2);

C2 = sqrt(g_GC * R_GC * T2);
Ct_2 = sqrt(C2^2 - Cz_2^2); % Velocità tangenziale

alpha2 = atan(Ct_2 / Cz_2);
alpha2_deg = rad2deg(alpha2);

A2 = m_0 / (rho2 * C2);
r_t=0.487;
r_h_in = sqrt(r_t^2-A2/pi);
r_m_in=0.5*(r_t+r_h_in);

%% Rotore

omega = 755.0118;   % da compressore

% Sezione d'ingresso
A2 = m_0 / (rho2 * Cz_2);
r_h2 = sqrt(r_t^2-A2/pi);
r_m2=(r_h2+r_t)*0.5;
U_m = omega * r_m;

W2_z = C2_z;
W2_t = Ct_2 - U_m;
W2 = sqrt(W2_t^2 + W2_z^2);
a2 = sqrt(g_GC * R_GC * T2);
% Fino a qua giusto
Tt2 = T_t4;
pt2 = P_t4;

% Sezione d'uscita
M_w3 = 0.9;
W3 = sqrt((a2^2 + (g_GC - 1) * (W2^2 / 2)) / ((1 / M_w3^2) + (g_GC - 1) / 2));

W3_z = C2_z;
W3_theta = sqrt(W3^2 - W3_z^2);

C3_t = U_m - W3_theta;
C3 = sqrt(C3_t^2 + W3_z^2);

R = 1 - (Ct_2 + C3_t) / (2 * U_m); % Grado di reazione

w_stadio = U_m * (Ct_2 - C3_t);
T_t3 = T_t2 - (w_stadio / cp_GC);

% Temperatura statica alla sezione 3
T3 = T_t3 - (C3^2) / (2 * cp_GC);

% Numero di Mach in uscita
a3 = sqrt(g_GC * R_GC * T3);
M3 = C3 / a3;

% Pressioni alla sezione 3
p_t3 = p_t2 * (T_t3 / T_t2)^(g_GC / ((g_GC - 1) * e_t));
p3 = p_t3 / (1 + (g_GC - 1)/2 * M3^2)^(g_GC / (g_GC - 1));

% Densità e Area di passaggio alla sezione 3
rho3 = p3 / (R_GC * T3);
A3 = m_g / (rho3 * W3_z);

% Geometria della pala in uscita
h3 = A3 / (2 * pi * r_m);
r_t3 = r_m + h3 / 2;
r_h3 = r_m - h3 / 2;

%% Sezione di scarico

T_t3=
T_t2=
Pc = ;
Pt = Pc / eta_m;

T_t5cool = ((1 - eps_cool) * cp_GC * T_t4 + eps_cool * cp_a * T_t3 - Pt/ m0) / ((1 - eps_cool) * cp_GC + eps_cool * cp_a);
taut = T_t5cool / T_t4;
pi_t = taut ^ (g_GC / ((g_GC - 1) * e_t));

T5_cool = T_t5cool - (0.5 * Cz_2^2) / cp_GC;
M5 = Cz_2 / sqrt(g_GC * R_GC * T5_cool);
P_t5 = pi_t * P_t4;
P5 = P_t5 / (1 + (M5^2) * (g_GC - 1) / 2);
rho5 = P5 / (R_GC * T5_cool);
A5 = m0 / (rho5 * Cz_2);

r_h_exit = 0.5*(r_t^2-A5/pi);
h5 = r_t-r_h_exit;
