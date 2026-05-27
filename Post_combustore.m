function [pf, Tf, Mf, A_ram] = Post_combustore()

%% Dimensionamento post combustore

supCruise = Supercruise_parameters();
Air = Air_parameters('25000');

m_a = supCruise.m_a;
f = supCruise.f;

% Parametri termodinamici
g_a = Air.g;
cp_a = Air.cp;              % Calore specifico aria [J/kgK]
cp_GC = Air.cp_GC;          % Calore specifico gas combusti [J/kgK]
H_f = 43000000;             % Potere calorifico [J/kg]
eta_AB = 0.98;              % Rendimento post-combustore
pi_AB = 0.97;

% Ingresso
[p_i, T_i, rho_i, M_i, A_i] = Presa('sup');

Ttot1 = T_i*(1+M_i^2*(g_a-1)/2);
ptot1 = p_i*(1+M_i^2*(g_a-1)/2)^(g_a/(g_a-1));
r_t5 = 0.492;     % Raggio tip in turbina (dovrebbe essere giusto)
L_H_ratio = 3;        %arbitrario
dhratio = 0.3;        %arbitrario

% Temperatura totale all'uscita
T_t7 = (cp_a * Ttot1 + f * H_f * eta_AB) / ((1 + f) * cp_GC);
p_t7 = ptot1*pi_AB;
Tf = T_t7/(1+M_i^2*(g_a-1)/2);
pf = p_t7 / (1+M_i^2*(g_a-1)/2)^(g_a/(g_a-1));
Mf = M_i;

% dimensioni post bruciatore
r_ram = sqrt((A_i / pi) + r_t5^2);
A_ram = p_i * r_ram^2;
H_ram = 2 * r_ram;
L_ram = L_H_ratio * H_ram;

% dimensioni stabilizzatori di fiamma (ne metto 8)
d_fh = dhratio * H_ram / 8;

