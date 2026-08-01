%% Dati

Air = Air_parameters('12000');
subCruise = Subcruise_parameters();

m_a = subCruise.m_a;
f = subCruise.f;
m_tot = m_a * (1 + f);
v0 = subCruise.v0;

p0 = Air.p;
g_GC = Air.g_GC;
R_GC = Air.R_GC;
cp_GC = Air.cp_GC;

[p_i, T_i, M_i, A_i] = Turbina();
p_t7 = p_i*(1+M_i^2*(g_GC-1)/2)^(g_GC/(g_GC-1));
T_t7 = T_i*(1+M_i^2*(g_GC-1)/2);
H_ram = 2*sqrt(A_i/pi);

%% Ugello supersonico

% Area di gola
% Funzione di Vandenkerckhove
Gamma = sqrt(g_GC) * (2 / (g_GC + 1))^((g_GC + 1) / (2 * (g_GC - 1)));
% Dimensionamento della gola
A_g = (m_tot * sqrt(R_GC * T_t7)) / (p_t7 * Gamma);

%% Ugello di de Laval supersonico

% Espansione isentropica da condizioni TOTALI
Te   = T_t7 * (p0 / p_t7)^((g_GC-1)/g_GC);
v_e  = sqrt(2 * cp_GC * (T_t7 - Te));

M_e  = v_e / sqrt(g_GC * R_GC * Te);

% Spinta netta (corretta se p_e = p0)
Effective_thrust = m_tot * v_e - m_a * v0;

r_g=sqrt(A_g/pi);
ratio_A = (1 / M_e) * ( (2 / (g_GC + 1)) * (1 + ((g_GC - 1) / 2) * M_e^2) )...
    ^((g_GC + 1) / (2 * (g_GC - 1)));
A_e = A_g * ratio_A;
r_e=sqrt(A_e/pi);
a1=deg2rad(40);
a2=deg2rad(15);
r_AB = 0.5763;
Lcon= (r_AB-r_g)/tan(a1);
Ldiv = (r_e-r_g)/tan(a2);
Ln=Lcon+Ldiv;