clear
close
clc

%% Dati

Air = Air_parameters('25000');
supCruise = Supercruise_parameters();

m_a = supCruise.m_a;
f = supCruise.f;
m_tot = m_a * (1 + f);

p0 = Air.p;
g_GC = Air.g_GC;
R_GC = Air.R_GC;
cp_GC = Air.cp_GC;

p_t7 = 1.3483e+05;
T_t7 = 1.461e+03;
H_ram = 1.6902;

soglia_n = 1.05;

%% Ugello supersonico

% Area di gola
% Funzione di Vandenkerckhove
Gamma = sqrt(g_GC) * (2 / (g_GC + 1))^((g_GC + 1) / (2 * (g_GC - 1)));

% Dimensionamento della gola
A_g = (m_tot * sqrt(R_GC * T_t7)) / (p_t7 * Gamma);

% Spinta ugello CD (espansione ottima)
v_e = sqrt(2 * cp_GC * T_t7 * (1 - (p0 / p_t7)^((g_GC-1)/g_GC)));
F_CD = m_tot * v_e;

% Spinta ugello Convergente (blocco sonico in gola)
% Condizioni critiche in gola
T_g = T_t7 * (2 / (g_GC + 1));
p_g = p_t7 * (2 / (g_GC + 1))^(g_GC / (g_GC - 1));
v_g = sqrt(g_GC * R_GC * T_g);
F_C = (m_tot * v_g) + A_g * (p_g - p0);

% Rapporto di Spinta al punto operativo
Ratio_F = F_CD / F_C;

if Ratio_F > soglia_n
    fprintf('Vantaggio del %.1f%% (> %.1f%%).\n', (Ratio_F-1)*100, (soglia_n-1)*100);
    disp('Ugello Convergente-Divergente richiesto.');
else
    fprintf('Vantaggio del %.1f%% (< %.1f%%).\n', (Ratio_F-1)*100, (soglia_n-1)*100);
    disp('Ugello Convergente sufficiente.');
end

%% Ugello di de Laval supersonico

% Mach di uscita
M_e = sqrt( (2 / (g_GC - 1)) * ( (p_t7 / p0)^((g_GC - 1) / g_GC) - 1 ) );

% Area di uscita
ratio_A = (1 / M_e) * ( (2 / (g_GC + 1)) * (1 + ((g_GC - 1) / 2) * M_e^2) )...
    ^((g_GC + 1) / (2 * (g_GC - 1)));
A_e = A_g * ratio_A;
D_g = 2 * sqrt(A_g / pi);
D_e = 2 * sqrt(A_e / pi);

% Scelta degli angoli
theta_c_deg = 50;
theta_d_deg = 25;

% Conversione in radianti
theta_c = deg2rad(theta_c_deg);
theta_d = deg2rad(theta_d_deg);

% Calcolo delle lunghezze geometriche
L_c = (H_ram - D_g) / (2 * tan(theta_c)); % Lunghezza del convergente
L_d = (D_e - D_g) / (2 * tan(theta_d));   % Lunghezza del divergente
L_tot = L_c + L_d;                        % Lunghezza totale dell'ugello
