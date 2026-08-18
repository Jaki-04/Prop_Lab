function [ptot_Out_cc, Ttot_Out_cc, MOut] = Camera()


S = Subcruise_parameters();
m_a = S.m_a;
options = optimoptions('fsolve','Display','none');
[ptot_In, Ttot_In, MIn, AIn] = Compressore();

% Diffusore prima della camera

G = Geometry();
rIn = G.Rc_he;
RIn = G.Rc_t;

Air = Air_parameters();
cp_a = Air.cp;
g_a = Air.g;
R_a = Air.R;
T_In = Ttot_In/(1+MIn^2*(g_a-1)/2);
p_In = ptot_In/(1+MIn^2*(g_a-1)/2)^(g_a/(g_a-1));


M_Out_diff = 0.2;             % Mach target
pi_diff = 0.97;
ptot_Out_diff = ptot_In * pi_diff;
p_Out_diff = ptot_Out_diff/(1+M_Out_diff^2*(g_a-1)/2)^(g_a/(g_a-1));
Ttot_Out_diff = Ttot_In;
T_Out_diff = Ttot_Out_diff/(1+M_Out_diff^2*(g_a-1)/2);
rho_Out_diff = p_Out_diff/(R_a*T_Out_diff);

AOut_diff = m_a/(rho_Out_diff*M_Out_diff*sqrt(g_a*R_a*T_Out_diff));
deltar = fsolve(@(dr) pi*((RIn)^2-(rIn-dr)^2)-AOut_diff, 0, options);
ROut_diff = RIn;
rOut_diff = rIn-deltar;
teta = deg2rad(3.5);
L_diff_cc = deltar/sin(teta);
WLratio_diff = L_diff_cc/(RIn-rIn);

% Camera di combustione
eps_cool = 0.0245;            % Spillamento di progetto
H_f = 43*10^6;
eta_b = 0.98;
MOut = 0.5;     % Mach in ingresso alla turbina

% Temperatura totale in uscita
g_GC = Air.g_GC;
cp_GC = Air.cp_GC;
Ttot_Out_cc =@(f) ((1-eps_cool)*cp_a*Ttot_Out_diff + f * H_f*eta_b)./((1+f-eps_cool)*cp_GC);
f=fsolve(@(f) Ttot_Out_cc(f) - 1420, 0.0294);
Ttot_Out_cc = Ttot_Out_cc(f);
T_Out_cc = Ttot_Out_cc/(1+MOut^2*(g_a-1)/2);

% Valori (dalle slide)
DP_P = 0.06;
DP_qref = 20;
I = 10*10^6;        % Capire che valore mettere (cercare qualche paper, slide dicono da 8 a 10)

Aref = ((R_a/2) * (m_a * T_In^0.5 / p_Out_diff)^2 * DP_qref * DP_P^-1)^0.5;
rm = (ROut_diff+rOut_diff)/2;
href = Aref/(2*pi*rm);
n = 1.8;                                                           % Esponente per l'efficienza cinetica (citare qualche paper)
Vcc = (H_f*f*m_a*(1-eps_cool)*eta_b)/(I*(p_Out_diff/10^5)^n);    % Equazione di Lefebvre per l'intensità di combustione
Lcc = Vcc/(2*pi*rm*href);                                       % Stima della lunghezza della camera

% Perdite di pressione Camera 
DP_qref = 20;
DP_P = 0.06;

rho_In = p_In / (R_a * T_In); % Densità ingresso diffusore
 
U_ref = m_a/(rho_In*Aref); % Velocità di riferimento 
q_ref = 0.5*rho_In*U_ref^2; % QdM di riferimento 
M_ref = U_ref/(g_a*R_a*T_In)^0.5; % Mach di riferimento

Dp_cold = 0.06 * ptot_In;
Dp_hot = q_ref * (T_Out_cc/T_In-1);
ptot_Out_cc = ptot_In - Dp_cold - Dp_hot;
p_Out_cc = ptot_Out_cc/(1+MOut^2*(g_a-1)/2)^(g_a/(g_a-1));

% Dp_cold_verifica = DP_qref+0.5*R_a*((m_a*T_In^0.5)/(Aref*p_In))^2*p_In;
