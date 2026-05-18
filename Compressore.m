%% Dimensionamento del compressore  (ford Kompressor)

%function [p_tOut, T_tOut, M_Out, A_Out_Comp, L_compressore] = Compressore()
Air = Air_parameters();
R = Air.R;
g = Air.g;
cp = Air.cp;
options = optimoptions('fsolve','Display','none');
S = Subcruise_parameters();
m_a = S.m_a;

[pf, Tf, rhof, Mf, Af] = Presa('sub');
% Valori di ingresso
    M_In = Mf;              
    p_In = pf;
    p_tIn = p_In*(1 + (g-1)/2 * M_In^2)^(g/(g-1));
    T_In = Tf;

% Determinazione sezione di uscita

C_z1 = M_In*sqrt(g*R*T_In);                              % Velocità assiale (costante)
A_in_comp = Af;                                         % Area di ingresso nel compressore

eta_poli=0.9;                                           % Rendimento politropico del compressore (costante)
beta = S.b;                                             % Rapporto di compressione (noto)

% Dimensionamento con primo stadio r_h1/r_t = 0.5
a=0.5;
r_t = sqrt(A_in_comp/(pi*(1-a^2)));
r_h1 = r_t*a;
r_pitch1 = (r_t+r_h1)/2;

% Omega necessarie ad avere Mrel=1 sulla pitchline del primo stadio
U_pitch = @(omega, r_pitch) omega*r_pitch;
Mrel = @(omega, r_pitch, T) sqrt(C_z1^2+U_pitch(omega, r_pitch)^2)/sqrt(g*R*T);
omega = fsolve(@(omega) Mrel(omega, r_pitch1, T_In)-1, 7000, options);

% Ciclo su tutti gli stadi necessari
A_In_stadio = A_in_comp;
T_In_stadio = T_In;
p_tIn_stadio = p_tIn;
M_In_stadio = M_In;
beta_real = 1;
N_stadi = 0;

l_stadio = zeros(1, 10);
stage_space = zeros(1, 10);
r_hub = zeros(1, 10);
cord_r = zeros(1, 10);
spacing_r = zeros(1, 10);
cord_s = zeros(1, 10);
spacing_s = zeros(1, 10);
bladeN_r = zeros(1, 10);
bladeN_s = zeros(1, 10);
gR = zeros(1, 10);
Dr = zeros(1, 10);
Ds = zeros(1, 10);

gamma_vec = zeros(10, 1);

while beta_real<beta
    N_stadi = N_stadi+1;
    r_h_stadio = sqrt(r_t^2-A_In_stadio/pi);
    r_hub(N_stadi) = r_h_stadio; 
    
    % Triangolo di velocità per lo stadio corrente
    r_pitch_stadio = (r_t+r_h_stadio)/2;
    U_pitch_stadio = omega*r_pitch_stadio;
    W1 = sqrt(U_pitch_stadio.^2 + C_z1.^2);
    W2 = 0.72*W1;
    gamma = asin(C_z1/W2);
    gamma_vec(N_stadi) = pi/2 - gamma;
    C2 = sqrt(W2^2+U_pitch_stadio^2-2*U_pitch_stadio*W2*cos(gamma));
    Cteta2 = sqrt(C2^2-C_z1^2);
    
    % Lavoro dello stadio corrente 
    L_stadio = U_pitch_stadio*(Cteta2);

    % Temperatura totale in uscita
    T_tIn_stadio = T_In_stadio*(1+M_In_stadio^2*(g-1)/2);
    T_Out_stadio = L_stadio/cp + T_In_stadio;
    M_Out_stadio = C_z1 / sqrt(g * R * T_Out_stadio); 
    T_tOut_stadio = T_Out_stadio*(1+M_Out_stadio^2*(g-1)/2);

    % Rapporto di compressione ottenuto (con rendimento)
    b_stadio = (T_tOut_stadio/T_tIn_stadio)^( 1/((g-1)/(eta_poli*g)) );
    p_tOut_stadio = b_stadio * p_tIn_stadio;                                                 
    p_Out_stadio = p_tOut_stadio/(1 + (g-1)/2 * M_Out_stadio^2)^(g/(g-1)); 

    % Sezione di uscita dello stadio
    rho_Out_stadio = p_Out_stadio/(R * T_Out_stadio);                            
    A_Out_stadio = m_a/(rho_Out_stadio * C_z1);                      
    
    
    % Coefficienti di diffusione (devono essere <=0.6 per non avere separazione)
    %sigma_r = 1;
    %sigma_s = 1.25;
    
    Dcoeff_r_stadio = @(sigma_r) 1-W2/W1+Cteta2/(2*sigma_r*W1);
    sigma_r = fsolve(@(sigma_r) Dcoeff_r_stadio(sigma_r) -0.55, 1, options);      % Impongo D=0.55 e trovo sigma_r
    Dr(N_stadi) = sigma_r;
    Dcoeff_s_stadio = @(sigma_s) 1-C_z1/C2+Cteta2/(2*sigma_s*C2);
    sigma_s = fsolve(@(sigma_s) Dcoeff_s_stadio(sigma_s) -0.55, 1, options);      % Impongo D=0.55 e trovo sigma_s
    Ds(N_stadi) = sigma_s;
    
    % Grado di reazione
    R_c_stadio = 1-Cteta2/(2*U_pitch_stadio);
    gR(N_stadi) = R_c_stadio;
    
    % Corda rotore e statore fissando il Re minimo
    nu1 = 4.67*10^-5;
    
    cr_stadio = 300000*nu1/W1;
    cord_r(N_stadi) = cr_stadio;
    
    cs_stadio = 300000*nu1/C2;
    cord_s(N_stadi) = cs_stadio;
    
    % Spaziatura pale rotore e statore usando i sigma dati
    
    sr_stadio = cr_stadio/sigma_r;
    spacing_r(N_stadi) = sr_stadio;
    
    ss_stadio = cs_stadio/sigma_s;
    spacing_s(N_stadi) = ss_stadio;

    % Lunghezza approssimativa dello stadio 
    l_r = cr_stadio*cos(pi/2-gamma);                                                % Corda inclinata uguale a W2;
    l_s = cs_stadio/2*sin(atan(C_z1/(U_pitch_stadio-W2*cos(gamma))))+cs_stadio/2;   % Corda metà inclinata come C2 e metà orizzontale
    stage_space(N_stadi) = cs_stadio*0.35;      % https://journals.sagepub.com/doi/10.1177/0957650914531949
    l_stadio(N_stadi) = (l_s+l_r+stage_space(N_stadi));        
     
    % Numero di pale
    
    Nr_stadio = (2*pi*r_pitch_stadio)/sr_stadio;
    bladeN_r(N_stadi) = ceil(Nr_stadio);
    
    Ns_stadio = (2*pi*r_pitch_stadio)/ss_stadio;
    bladeN_s(N_stadi) = ceil(Ns_stadio);
    
    % Aggiorno beta_real
    beta_real = beta_real*b_stadio;

    % Aggiorno quantità in ingresso allo stadio successivo
    A_In_stadio = A_Out_stadio;
    T_In_stadio = T_Out_stadio;
    p_tIn_stadio = p_tOut_stadio;
    M_In_stadio = M_Out_stadio;
    
end

%svergolamento palette
Inc = atan((omega*r_pitch1)/C_z1) - gamma_vec(1);
gamma_sverg = @(r) Inc - atan((omega*r)/C_z1);


p_tOut = p_tIn_stadio;
M_Out = M_In_stadio;
T_tOut = T_In_stadio*(1+M_Out^2*(g-1)/2);
A_Out_Comp = A_In_stadio;
L_compressore= sum(l_stadio)+stage_space(end);