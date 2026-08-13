function [p_t7, T_t7, Mf, A_AB] = Post_combustore()

%% Dimensionamento post combustore

supCruise = Supercruise_parameters();
Air = Air_parameters('25000');

m_a = supCruise.m_a;
f = supCruise.f;

% Parametri termodinamici
g_a = Air.g;
g_GC = Air.g_GC;
cp_a = Air.cp;              % Calore specifico aria [J/kgK]
R_a = Air.R; 
cp_GC = Air.cp_GC;          % Calore specifico gas combusti [J/kgK]
R_GC = Air.R_GC;
H_f = 43000000;             % Potere calorifico [J/kg]
eta_AB = 0.98;              % Rendimento post-combustore
pi_AB = 0.97;

% Ingresso
[p_i, T_i, rho_i, M_i, A_i] = Presa('sup');

Ttot1 = T_i*(1+M_i^2*(g_a-1)/2);
ptot1 = p_i*(1+M_i^2*(g_a-1)/2)^(g_a/(g_a-1));
r_t5 = 0.484;        % Raggio tip in turbina
L_H_ratio = 2.2;        %arbitrario
dhratio = 0.25;        %arbitrario

% Temperatura totale all'uscita
T_t7 = (cp_a * Ttot1 + f * H_f * eta_AB) / ((1 + f) * cp_GC);
p_t7 = ptot1*pi_AB;
Mf=0.469;


% dimensioni condotti di bypass
r_ram = sqrt((A_i / pi) + r_t5^2);

% Diffusore subsonico prima dell'AB
M2 = 0.25;                                  % Mach target all'ingresso del postbruciatore 
T2 = Ttot1/(1+M2^2*(g_a-1)/2);
v2 = M2*sqrt(R_a*g_a*T2);
p2 = ptot1/((1+M2^2*(g_a-1)/2)^(g_a/(g_a-1)));
rho2 = p2/(R_a*T2); 
A_AB = m_a/(rho2*v2); 

r_AB = sqrt(A_AB/pi);
H_AB = 2 * r_AB;
L_AB= L_H_ratio * H_AB;

% dimensioni stabilizzatori di fiamma (ne metto 2)
d_fh = dhratio * H_AB / (2*2);

