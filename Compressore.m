%% Dimensionamento del compressore

Air = Air_parameters();
R=Air.R;
g = Air.g;

S = Subcruise_parameters();
m_a = S.m_a;

Min = 0.3;              % Mach di ingresso
p_tIn = 3.084344e+04;
p_In = 2.8976481e+04;
T_tIn = 247.955925;
T_In = 243.571636;
rho_In = p_In/(R*T_In);

C_z1 = 93.85278;                                                % Velocità assiale (costante)
A_in_comp = m_a/(rho_In * C_z1);                                %area di ingresso nel compressore

eta_poli=0.9;
beta = subCruise.b;
p_tOut = beta * p_tIn;
T_tOut = beta^((g-1)/(eta_poli*g)) * T_tIn;                          %temperatura totale in uscita con rendimento politropico
T_Out = T_tOut - (g-1)/(2*g*R)*C_z1^2;                              %temperatura statica in uscita
M_Out = C_z1 / sqrt(g * R * T_Out);                                %mach in uscita
p_Out = p_tOut/(1 + (g-1)/2 * M_Out^2)^(g/(g-1));                  %pressione statica dell'aria in uscita
rho_Out = p_Out/(R * T_Out);                                        %rho dell'aria in uscita
A_Out_comp = m_a/(rho_Out * C_z1);                                     %area di uscita del compressore

%dimensionamento primo stadio r_h1/r_t1 = 0.5
a=0.5;
r_t = sqrt(A_in_comp/(pi*(1-a^2)));
r_h1 = r_t*a;
r_he = sqrt(r_t^2-(A_Out_comp)/pi);
U_tip = 400;
omega = U_tip/r_t;

U_pala = @(r) (omega * r);                                                %r può variare solo tra r_h e r_t
W = @(U_pala) sqrt(U_pala.^2 + (C_z1 .* ones(size(U_pala))).^2);
Re_corda = @(chord, W) (rho_in * chord .* W) / 10^(-3);        %fcn handle più generica per il reynolds sulla pala
%disp(Re_corda(0.2, U_pala(r_h1:0.01:r_t)));
disp(omega/(2*pi) * 60);
%disp(U_tip/sqrt(g*R*T2_sub))