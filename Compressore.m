%% Dimensionamento del compressore

Air = Air_parameters();
R = Air.R;
g = Air.g;
cp = Air.cp;

S = Subcruise_parameters();
m_a = S.m_a;

% Valori di ingresso
    Min = 0.5;              
    p_tIn = 3.049306402121894e+04;
    p_In = p_tIn/(1 + (g-1)/2 * Min^2)^(g/(g-1));
    T_tIn = 247.955925;
    T_In = T_tIn/(1 + (g-1)/2 * Min^2);
    rho_In = p_In/(R*T_In);

C_z1 = Min*sqrt(g*R*T_In);                              % Velocità assiale (costante)
A_in_comp = m_a/(rho_In * C_z1);                        % Area di ingresso nel compressore

eta_poli=0.9;                                           % Rendimento politropico del compressore (costante)
beta = S.b;                                             % Rapporto di compressione (noto)
lpbeta = 4;
beta1 = lpbeta;
beta2 = beta/lpbeta;

p_tmed = beta1 * p_tIn;                                  % Pressione totale in uscita
T_tmed = beta1^((g-1)/(eta_poli*g)) * T_tIn;             % Temperatura totale in uscita con rendimento politropico
T_med = T_tmed - (g-1)/(2*g*R)*C_z1^2;                  % Temperatura statica in uscita
M_med = C_z1 / sqrt(g * R * T_med);                     % Mach in uscita
p_med = p_tmed/(1 + (g-1)/2 * M_med^2)^(g/(g-1));       % Pressione statica dell'aria in uscita
rho_med = p_med/(R * T_med);                            % Rho dell'aria in uscita
A_med_comp = m_a/(rho_med * C_z1);                      % Area di uscita del compressore



p_tOut = beta2 * p_tmed;                                  % Pressione totale in uscita
T_tOut = beta2^((g-1)/(eta_poli*g)) * T_tmed;             % Temperatura totale in uscita con rendimento politropico
T_Out = T_tOut - (g-1)/(2*g*R)*C_z1^2;                  % Temperatura statica in uscita
M_Out = C_z1 / sqrt(g * R * T_Out);                     % Mach in uscita
p_Out = p_tOut/(1 + (g-1)/2 * M_Out^2)^(g/(g-1));       % Pressione statica dell'aria in uscita
rho_Out = p_Out/(R * T_Out);                            % Rho dell'aria in uscita
A_Out_comp = m_a/(rho_Out * C_z1);                      % Area di uscita del compressore

% Dimensionamento con primo stadio r_h1/r_t1 = 0.5
a=0.5;
r_t = sqrt(A_in_comp/(pi*(1-a^2)));
r_h1 = r_t*a;
r_HP = sqrt(r_t^2-(A_med_comp)/pi);
r_he = sqrt(r_t^2-(A_Out_comp)/pi);
r_pitch1 = (r_t+r_h1)/2;
r_pitch2 = (r_t+r_HP)/2;

% Omega necessarie ad avere Mrel=1 sulla pitchline
U_pitch = @(omega, r_pitch) omega*r_pitch;
Mrel = @(omega, r_pitch, T) sqrt(C_z1^2+U_pitch(omega, r_pitch)^2)/sqrt(g*R*T);
omega1 = fsolve(@(omega) Mrel(omega, r_pitch1, T_In)-1, 7000);
omega2 = fsolve(@(omega) Mrel(omega, r_pitch2, T_med)-1, 7000);


U_pitch1 = omega1*r_pitch1;
W1_LP = sqrt(U_pitch1.^2 + C_z1.^2);
W2_LP = 0.72*W1_LP;
gammaLP = asin(C_z1/W2_LP);
C2_LP = sqrt(W2_LP^2+U_pitch1^2-2*U_pitch1*W2_LP*cos(gammaLP));
Cteta2_LP = sqrt(C2_LP^2-C_z1^2);

U_pitch2 = omega2*r_pitch2;
W1_HP = sqrt(U_pitch2.^2 + C_z1.^2);
W2_HP = 0.72*W1_HP;
gammaHP = asin(C_z1/W2_HP);
C2_HP = sqrt(W2_HP^2+U_pitch2^2-2*U_pitch2*W2_HP*cos(gammaHP));
Cteta2_HP = sqrt(C2_HP^2-C_z1^2);

% Coefficienti di diffusione (devono essere <=0.6 per non avere separazione)
sigma_r = 1;
sigma_s = 1.25;

Dcoeff_LPr = 1-W2_LP/W1_LP+Cteta2_LP/(2*sigma_r*W1_LP);
Dcoeff_LPs = 1-C_z1/C2_LP+Cteta2_LP/(2*sigma_s*C2_LP);

Dcoeff_HPr = 1-W2_HP/W1_HP+Cteta2_HP/(2*sigma_r*W1_HP);
Dcoeff_HPs = 1-C_z1/C2_HP+Cteta2_HP/(2*sigma_s*C2_HP);

% Grado di reazione
R_LP = 1-Cteta2_LP/(2*U_pitch1);
R_HP = 1-Cteta2_HP/(2*U_pitch1);

% Corda rotore e statore
nu1 = 4.67*10^-5;

cr_LP = 300000*nu1/W1_LP;
cr_HP = 300000*nu1/W1_HP;

cs_LP = 300000*nu1/C2_LP;
cs_HP = 300000*nu1/C2_HP;

% Spaziatura rotore e statore

sr_LP = cr_LP/sigma_r;
sr_HP = cr_HP/sigma_r;

ss_LP = cs_LP/sigma_s;
ss_HP = cs_HP/sigma_s;

% Numero di pale

Nr_LP = (2*pi*r_pitch1)/sr_LP;
Nr_HP = (2*pi*r_pitch2)/sr_HP;

Ns_LP = (2*pi*r_pitch1)/ss_LP;
Ns_HP = (2*pi*r_pitch2)/ss_HP;

% Lavoro di uno stadio

L_LP = U_pitch1*(Cteta2_LP);
Ltot_LP = cp*(T_med-T_In);
N_stadiLP = Ltot_LP/L_LP;
b_stadio_LP = beta1^(1/N_stadiLP);
N_stadiLP = ceil(N_stadiLP);
beta1_real = b_stadio_LP^N_stadiLP;

L_HP = U_pitch2*(Cteta2_HP);
Ltot_HP = cp*(T_Out-T_med);
N_stadiHP = Ltot_HP/L_HP;
b_stadio_HP = beta2^(1/N_stadiHP);
N_stadiHP = ceil(N_stadiHP);
beta2_real = b_stadio_HP^N_stadiHP;

beta = beta1_real*beta2_real;