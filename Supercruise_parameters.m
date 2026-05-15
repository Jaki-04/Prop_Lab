%% Studio parametrico - Crociera supersonica

function supCruise = Supercruise_parameters(varargin)

% Rendimenti e parametri
pi_AB = 0.95;       % Fissato
pi_d_sup = 0.74;     % Fissato (presa obiettivo)
eta_AB=0.92;
eta_n = 0.98;
H_f=43000000;       % Fissato
Tmax_AB = 2100;

% Aria a 25km
Air = Air_parameters('25000');
p=Air.p;
T=Air.T;
cp_a=Air.cp;
g_a = Air.g;
R_a = Air.R;


% Caratteristiche GC
cp_GC=Air.cp_GC;
g_GC =Air.g_GC;

% Ciclo su Ttot_AB
Ttot_AB = 

% Crociera supersonica
M_supercruise = 3.5;
v0_supercruise = M_supercruise*sqrt(g_a*R_a*T);

% Presa
Ttot1=T*(1+(g_a-1)/2*M_supercruise^2);
ptot1=p*(1+(g_a-1)/2*M_supercruise^2)^(g_a/(g_a-1))*pi_d_sup;

% AB
ptot2=pi_AB*ptot1;
f = (cp_a.*Ttot1-Ttot2.*cp_GC)./(Ttot2.*cp_GC-H_f.*eta_b);

% Ugello
T_ratio = (p/ptot2).^( (g_GC-1)/(g_GC) );
ve = sqrt( 2.*cp_GC.*Ttot2.*eta_n.*(1-T_ratio));

% Parametri di merito
I_sp_a = ((1+f).*ve-v0_supercruise);
TSFC=f./I_sp_a;

[TSFC_min, F_min] = min(TSFC);

if ismember('plot', varargin)
    figure()
    plot(f, TSFC*3600);
    hold on;
    yline(TSFC_min*3600, '--k');
    plot(f(F_min), TSFC_min*3600, 'o', 'MarkerFaceColor', 'r');
    ylim([0.001, 0.5])
    legend("$TSFC$","", "$TSFC_{min}$", 'Interpreter','latex');
end

supCruise.m_a=m_a(F_min);
supCruise.f = f(F_min);
supCruise.I_sp_a=I_sp_a(F_min);
supCruise.TSFC=TSFC(F_min);
supCruise.M=M_supercruise;
supCruise.v0=v0_supercruise;