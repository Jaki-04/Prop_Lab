clear
close
clc

%% Dati

Air = Air_parameters('25000');
supCruise = Supercruise_parameters();

m_a = supCruise.m_a;
f = supCruise.f;
m_tot = m_a * (1 + f);
v0 = supCruise.v0;

p0 = Air.p;
g_GC = Air.g_GC;
R_GC = Air.R_GC;
cp_GC = Air.cp_GC;

[p_ti, T_ti, M_i, A_i] = Post_combustore();
H_ram = 2*sqrt(A_i/pi);

%% Ugello supersonico

% Area di gola
% Funzione di Vandenkerckhove
Gamma = sqrt(g_GC) * (2 / (g_GC + 1))^((g_GC + 1) / (2 * (g_GC - 1)));

% Dimensionamento della gola
A_g = (m_tot * sqrt(R_GC * T_ti)) / (p_ti * Gamma);

%% Ugello di de Laval supersonico

% Espansione isentropica da condizioni TOTALI
p_t9=p_ti;
T_t9=T_ti;
T9  = @(p) T_t9 * (p / p_t9)^((g_GC-1)/g_GC);
v9  = @(p) sqrt(2 * cp_GC * (T_t9 - T9(p)));
M9  = @(p) v9(p) / sqrt(g_GC * R_GC * T9(p));
ratio_A = @(p) (1 / M9(p)) * ( (2 / (g_GC + 1)) * (1 + ((g_GC - 1) / 2) * M9(p)^2) )...
    ^((g_GC + 1) / (2 * (g_GC - 1)));
A_e = @(p) A_g * ratio_A(p);

% Spinta netta (corretta se p_e = p0)
Effective_thrust = @(p) m_tot * v9(p) - m_a * v0 + (p-p0)*A_e(p);

%pe=fsolve(@(p) Effective_thrust(p) - 64000, 5000);
pe = 2549;
Effective_thrust(p0)

% Area di uscita
T9 = T9(pe);
v9 = v9(pe);
M9 = M9(pe);
A_e = A_e(pe);

D_g = 2 * sqrt(A_g / pi);
D_e = 2 * sqrt(A_e / pi);
% Scelta degli angoli
theta_c_deg = 45;
theta_d_deg = 25;

% Conversione in radianti
theta_c = deg2rad(theta_c_deg);
theta_d = deg2rad(theta_d_deg);

% Calcolo delle lunghezze geometriche
L_c = (H_ram - D_g) / (2 * tan(theta_c)); % Lunghezza del convergente
L_d = (D_e - D_g) / (2 * tan(theta_d));   % Lunghezza del divergente
L_tot = L_c + L_d;                        % Lunghezza totale dell'ugello
