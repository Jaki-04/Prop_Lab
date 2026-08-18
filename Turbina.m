function [p_5, T_5cool, M5, A5] = Turbina()

subCruise = Subcruise_parameters();
Air = Air_parameters('12000');

m_a = subCruise.m_a;
f = 0.0285;
cp_GC = Air.cp_GC;
g_GC = Air.g_GC;
cp_a = Air.cp;
R_GC = Air.R_GC;

M4 = 0.5;
eta_m = 0.92;
eps_cool = 0.0245;
m_0 = m_a * (1 + f);
e_t = 0.9;

[P_t4, T_t4]=Camera();

%% Statore

% Ingresso statore
T4 = T_t4 / (1 + ((g_GC - 1) * (M4 ^ 2)) / 2);
P4 = P_t4 / ((1 + ((g_GC - 1) * (M4 ^ 2)) / 2)) ^ (g_GC / (g_GC - 1));
Ci = M4 * sqrt(g_GC * R_GC * T4);
rho4 = P4 / (R_GC * T4);
Ai = m_0 / (rho4 * Ci);

Mstar = 1.0;          
Cz_ii = Ci;         

% Termodinamica in uscita
Tii = T_t4* (2 / (g_GC + 1));
Pii = P_t4 * (2 / (g_GC + 1)) ^ (g_GC / (g_GC - 1));
rhoii = Pii / (R_GC * Tii);

Cii = sqrt(g_GC * R_GC * Tii);
Ct_ii = sqrt(Cii^2 - Cz_ii^2); % Velocità tangenziale

alpha_ii = atan(Ct_ii / Cz_ii);
alpha_ii_deg = rad2deg(alpha_ii);

Astar = m_0 / (rhoii * Cii);
r_t=0.4842;
r_h_in = sqrt(r_t^2-Ai/pi);
r_m_in=0.5*(r_t+r_h_in);

%% Rotore

omega = 735.2970;   % da compressore

% Sezione d'ingresso
Aii = m_0 / (rhoii * Cz_ii);
r_h2 = sqrt(r_t^2-Aii/pi);
r_m2=(r_h2+r_t)*0.5;
h2=r_t-r_h2;
U_m = omega * r_m2;

Wii_z = Cz_ii;
Wii_t = Ct_ii - U_m;
Wii = sqrt(Wii_t^2 + Wii_z^2);
aii = sqrt(g_GC * R_GC * Tii);

T_tii = T_t4;
p_tii = P_t4;

% Sezione d'uscita
M_wiii = 0.9;
Wiii = sqrt((aii^2 + (g_GC - 1) * (Wii^2 / 2)) / ((1 / M_wiii^2) + (g_GC - 1) / 2));

Wz_iii = Cz_ii;
Wt_iii = sqrt(Wiii^2 - Wz_iii^2);

Ct_iii = U_m - Wt_iii;
Ciii = sqrt(Ct_iii^2 + Wz_iii^2);
theta_iii=atan(Ct_iii/Wz_iii);
alpha_iii=atan(Wt_iii/Wz_iii);

R = 1 - (Ct_ii + Ct_iii) / (2 * U_m); % Grado di reazione

w_stadio = U_m * (Ct_ii - Ct_iii);
T_tiii = T_tii - (w_stadio / cp_GC);

% Temperatura statica alla sezione 3
Tiii = T_tiii - (Ciii^2) / (2 * cp_GC);

% Numero di Mach in uscita
a3 = sqrt(g_GC * R_GC * Tiii);
M3 = Ciii / a3;

% Pressioni alla sezione 3
p_tiii = p_tii * (T_tiii / T_tii)^(g_GC / ((g_GC - 1) * e_t));
piii = p_tiii / (1 + (g_GC - 1)/2 * M3^2)^(g_GC / (g_GC - 1));

% Densità e Area di passaggio alla sezione 3
rhoiii = piii / (R_GC * Tiii);
Aiii = m_0 / (rhoiii * Wz_iii);

% Geometria della pala in uscita
r_h3 = sqrt(r_t^2-Aiii/pi);
h3 = r_t-r_h3;
r_p3=(r_h3+r_t)*0.5;

%% Sezione di scarico

T_t3=581.1140;
T_t2=247.9559;
Pc = m_a*cp_a*(T_t3-T_t2);
Pt = Pc / eta_m;

T_t5cool = ((1 - eps_cool+f) * cp_GC * T_t4 + eps_cool * cp_a * T_t3 - Pt/ m_a) / ((1 - eps_cool+f) * cp_GC + eps_cool * cp_a);
taut = T_t5cool / T_t4;
pi_t = taut ^ (g_GC / ((g_GC - 1) * e_t));

T_5cool = T_t5cool - (0.5 * Cz_ii^2) / cp_GC;
M5 = Cz_ii / sqrt(g_GC * R_GC * T_5cool);
p_t5 = pi_t * P_t4;
p_5 = p_t5 / (1 + (M5^2) * (g_GC - 1) / 2)^(g_GC / (g_GC - 1));
rho5 = p_5 / (R_GC * T_5cool);
A5 = m_0 / (rho5 * Cz_ii);

r_h_exit = sqrt(r_t^2-A5/pi);

N_stadi=Pt/(w_stadio*m_0);

L_tb = 0.1181*2 + 1.5*(r_t-r_h_exit);