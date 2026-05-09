clc;
clearvars;


S = Subcruise_parameters();
m_a = S.m_a;
f = S.f;
options = optimoptions('fsolve','Display','none');
[ptot_In, Ttot_In, MIn, AIn] = Compressore();

%% Diffusore prima della camera

G = Geometry();
rIn = G.Rc_he;
RIn = G.Rc_t;

Air = Air_parameters();
cp_a = Air.cp;
g_a = Air.g;
R_a = Air.R;


M_Out_diff = 0.2;             % Mach target
pi_diff = 0.97;
ptot_Out_diff = ptot_In * pi_diff;
p_Out_diff = ptot_Out_diff/(1+M_Out_diff^2*(g_a-1)/2)^(g_a/(g_a-1));
Ttot_Out_diff = Ttot_In;
T_Out_diff = Ttot_Out_diff/(1+M_Out_diff^2*(g_a-1)/2);
rho_Out_diff = p_Out_diff/(R_a*T_Out_diff);

AOut_diff = m_a/(rho_Out_diff*M_Out_diff*sqrt(g_a*R_a*T_Out_diff));
deltar = fsolve(@(dr) pi*((RIn+dr)^2-(rIn-dr)^2)-AOut_diff, 0, options);
ROut_diff = RIn+deltar;
rOut_diff = rIn-deltar;
teta = deg2rad(4);
L_diff_cc = deltar/sin(teta);
WLratio_diff = L_diff_cc/(RIn-rIn);

%% Camera di combustione
eps_cool = 0.036;            % Spillamento di progetto
H_f = 43000000;
eta_b = 0.98;

% Temperatura in uscita
g_GC = Air.g_GC;
cp_GC = Air.g_GC;
T_Out_cc= ((1-eps_cool)*cp_a*T_Out_diff + f * H_f*eta_b)./((1+f-eps_cool)*cp_GC);

% Valori (dalle slide)
DP_P = 0.06;
DP_qref = 20;
I = 10*10^6;

Aref = ((R_a/2) * (m_a * T_Out_diff^0.5 / p_Out_diff)^2 * DP_qref * DP_P^-1)^0.5;
rm = (ROut_diff+rOut_diff)/2;
href = Aref/(2*pi*rm);
Vcc = (H_f*f*m_a*(1-eps))/(I*(p_Out_diff/10^5));
Lcc = Vcc/(2*pi*rm*href);
