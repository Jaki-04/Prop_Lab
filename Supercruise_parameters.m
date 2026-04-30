%% Studio parametrico - Crociera supersonica

function [m_a, f, I_sp_a, TSFC] = Supercruise_parameters(varargin)

% Rendimenti e parametri
pi_AB = 0.98;       % Fissato
pi_d_sup = 0.7;     % Fissato (presa obiettivo)
eta_AB=0.92;
eta_n = 0.92;
H_f=43000000;       % Fissato
Tmax_AB = 2100;

% Caratteristiche aria e GC
cp_GC=1150;
cp_a=1000;
gamma_a=1.4;
gamma_GC = 1.33;
R_a=287.01;

% Aria a 25km
p_25000= 2549; 
T_25000= 216.65;
rho_25000 = p_25000/(R_a*T_25000);

f=0.0045:0.0001:0.1;

% Crociera supersonica
M_supercruise = 3.5;
v0_supercruise = M_supercruise*sqrt(gamma_a*R_a*T_25000);
T_supercruise=65000;

% Presa
Ttot1=T_25000*(1+(gamma_a-1)/2*M_supercruise^2);
ptot1=p_25000*(1+(gamma_a-1)/2*M_supercruise^2)^(gamma_a/(gamma_a-1))*pi_d_sup;

% AB
ptot2=pi_AB*ptot1;
Ttot2=(cp_a*Ttot1+f*H_f*eta_AB)./((1+f)*cp_GC);

% Ugello
T_ratio = 1-(p_25000/ptot2).^( (gamma_GC-1)/(gamma_GC) );
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
    legend("$TSFC$", "$I_sp}$", "$ve$", 'Interpreter','latex');
end