%% Studio parametrico - Crociera supersonica

function supCruise = Supercruise_parameters(varargin)

% Rendimenti e parametri
pi_AB = 0.98;       % Fissato
pi_d_sup = 0.7;     % Fissato (presa obiettivo)
eta_AB=0.92;
eta_n = 0.92;
H_f=43000000;       % Fissato
Tmax_AB = 2100;

% Caratteristiche GC
cp_GC=1150;
g_GC = 1.33;


% Aria a 25km
Air = Air_parameters('25000');
p=Air.p;
T=Air.T;
cp_a=Air.cp;
g_a = Air.g;
R_a = Air.R;

% Ciclo su f
f=0.0045:0.0001:0.1;

% Crociera supersonica
M_supercruise = 3.5;
v0_supercruise = M_supercruise*sqrt(g_a*R_a*T);
T_supercruise=65000;

% Presa
Ttot1=T*(1+(g_a-1)/2*M_supercruise^2);
ptot1=p*(1+(g_a-1)/2*M_supercruise^2)^(g_a/(g_a-1))*pi_d_sup;

% AB
ptot2=pi_AB*ptot1;
Ttot2=(cp_a*Ttot1+f*H_f*eta_AB)./((1+f)*cp_GC);

% Ugello
T_ratio = 1-(p/ptot2).^( (g_GC-1)/(g_GC) );
ve = sqrt( 2.*cp_GC.*Ttot2.*eta_n.*T_ratio);

% Portata d'aria target e portata massica corrispondente
m_a=T_supercruise./((1+f).*ve-v0_supercruise);
m_f =f.*m_a;

% Parametri di merito

I_sp_a = T_supercruise./m_a;
TSFC=m_f./T_supercruise;
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