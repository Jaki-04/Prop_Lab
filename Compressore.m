%% Dimensionamento del compressore  (ford Kompressor)

function [p_tOut, T_tOut, M_Out, A_Out_Comp] = Compressore()
Air = Air_parameters();
R = Air.R;
g = Air.g;
cp = Air.cp;
options = optimoptions('fsolve','Display','none');
S = Subcruise_parameters();
m_a = S.m_a;

[pf, Tf, rhof, Mf, Af] = Presa('sub');
% Valori di ingresso
    Min = Mf;              
    p_In = pf;
    p_tIn = p_In*(1 + (g-1)/2 * Min^2)^(g/(g-1));
    T_In = Tf;
    T_tIn = T_In*(1 + (g-1)/2 * Min^2);
    rho_In = rhof;

C_z1 = Min*sqrt(g*R*T_In);                              % Velocità assiale (costante)
A_in_comp = Af;                                         % Area di ingresso nel compressore

eta_poli=0.9;                                           % Rendimento politropico del compressore (costante)
beta = S.b;                                             % Rapporto di compressione (noto)

p_tOut = beta * p_tIn;                                % Pressione totale in uscita
T_tOut = beta^((g-1)/(eta_poli*g)) * T_tIn;           % Temperatura totale in uscita con rendimento politropico
T_Out = T_tOut - (g-1)/(2*g*R)*C_z1^2;                  % Temperatura statica in uscita
M_Out = C_z1 / sqrt(g * R * T_Out);                     % Mach in uscita
p_Out = p_tOut/(1 + (g-1)/2 * M_Out^2)^(g/(g-1));       % Pressione statica dell'aria in uscita
rho_Out = p_Out/(R * T_Out);                            % Rho dell'aria in uscita
A_Out_comp = m_a/(rho_Out * C_z1);                      % Area di uscita del compressore HP

% Dimensionamento con primo stadio r_h1/r_t1 = 0.5
a=0.5;
r_t = sqrt(A_in_comp/(pi*(1-a^2)));
r_h1 = r_t*a;
r_he = sqrt(r_t^2-(A_Out_comp)/pi);
r_pitch = (r_t+r_h1)/2;

% Omega necessarie ad avere Mrel=1 sulla pitchline
U_pitch = @(omega, r_pitch) omega*r_pitch;
Mrel = @(omega, r_pitch, T) sqrt(C_z1^2+U_pitch(omega, r_pitch)^2)/sqrt(g*R*T);
omega = fsolve(@(omega) Mrel(omega, r_pitch, T_In)-1, 7000, options);

U_pitch = omega*r_pitch;
W1 = sqrt(U_pitch.^2 + C_z1.^2);
W2 = 0.72*W1;
gamma = asin(C_z1/W2);
C2 = sqrt(W2^2+U_pitch^2-2*U_pitch*W2*cos(gamma));
Cteta2 = sqrt(C2^2-C_z1^2);

% Coefficienti di diffusione (devono essere <=0.6 per non avere separazione)
sigma_r = 1;
sigma_s = 1.25;

Dcoeff_r = 1-W2/W1+Cteta2/(2*sigma_r*W1);
Dcoeff_s = 1-C_z1/C2+Cteta2/(2*sigma_s*C2);

% Grado di reazione
R_c = 1-Cteta2/(2*U_pitch);

% Corda rotore e statore
nu1 = 4.67*10^-5;

cr = 300000*nu1/W1;

cs = 300000*nu1/C2;

% Spaziatura rotore e statore

sr = cr/sigma_r;

ss = cs/sigma_s;

% Numero di pale

Nr = (2*pi*r_pitch)/sr;

Ns = (2*pi*r_pitch)/ss;

% Lavoro di uno stadio

L_comp = U_pitch*(Cteta2);
Ltot_comp = cp*(T_Out-T_In);
N_stadi = Ltot_comp/L_comp;
b_stadio = beta^(1/N_stadi);
N_stadi = ceil(N_stadi);
beta_real = b_stadio^(N_stadi);

%calcolo le quantità con il beta reale (più realistico wow) 

p_tOut = beta_real * p_tIn;                                % Pressione totale in uscita
T_tOut = beta_real^((g-1)/(eta_poli*g)) * T_tIn;           % Temperatura totale in uscita con rendimento politropico
T_Out = T_tOut - (g-1)/(2*g*R)*C_z1^2;                  % Temperatura statica in uscita
M_Out = C_z1 / sqrt(g * R * T_Out);                     % Mach in uscita
p_Out = p_tOut/(1 + (g-1)/2 * M_Out^2)^(g/(g-1));       % Pressione statica dell'aria in uscita
rho_Out = p_Out/(R * T_Out);                            % Rho dell'aria in uscita
A_Out_Comp = m_a/(rho_Out * C_z1);                      % Area di uscita del compressore HP
