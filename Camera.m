clc;
clearvars;


S = Subcruise_parameters();
m_a = S.m_a;
options = optimoptions('fsolve','Display','none');
[ptot_In, Ttot_In, MIn, AIn] = Compressore();

%% Diffusore prima della camera

G = Geometry();
rIn = G.Rc_he;
RIn = G.Rc_t;

Air = Air_parameters();
g = Air.g;
R = Air.R;

M_Out = 0.2;             % Mach target
pi_diff = 0.97;
ptot_Out = ptot_In * pi_diff;
p_Out = ptot_Out/(1+M_Out^2*(g-1)/2)^(g/(g-1));
Ttot_Out = Ttot_In;
T_Out = Ttot_Out/(1+M_Out^2*(g-1)/2);
rho_Out = p_Out/(R*T_Out);

AOut = m_a/(rho_Out*M_Out*sqrt(g*R*T_Out));
deltar = fsolve(@(dr) pi*((RIn+dr)^2-(rIn-dr)^2)-AOut, 0, options);
ROut = RIn+deltar;
rOut = rIn-deltar;
teta = deg2rad(3);
L_diff_cc = deltar/sin(teta);

%% Camera di combustione
