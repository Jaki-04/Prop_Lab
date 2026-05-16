%% Studio parametrico - Crociera supersonica

function supCruise = Supercruise_parameters(varargin)

% Rendimenti e parametri
pi_AB = 0.95;       % Fissato
pi_d_sup = 0.73;     % Fissato (presa obiettivo)
eta_AB=0.98;
eta_n = 0.98;
H_f=43000000;       % Fissato
Tmax_AB = 2100;
T_supercruise = 65000;

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
R_GC = Air.R_GC;

% Ciclo su Ttot_AB
Ttot2 = [700:5:790, 791:0.1:801, 802:1:Tmax_AB];

% Crociera supersonica
M_supercruise = 3.5;
v0_supercruise = M_supercruise*sqrt(g_a*R_a*T);

% Presa
Ttot1=T*(1+(g_a-1)/2*M_supercruise^2);
ptot1=p*(1+(g_a-1)/2*M_supercruise^2)^(g_a/(g_a-1))*pi_d_sup;

% AB
ptot2=pi_AB*ptot1;
f = (cp_a.*Ttot1-Ttot2.*cp_GC)./(Ttot2.*cp_GC-H_f.*eta_AB);

% Ugello
T_ratio = (p/ptot2).^( (g_GC-1)/(g_GC) );
ve = sqrt( 2.*cp_GC.*Ttot2.*eta_n.*(1-T_ratio));
pe = p;
Te = T_ratio.*Ttot2;
rhoe = pe./(R_GC*Te);

% Parametri di merito
I_sp_a = ((1+f).*ve-v0_supercruise);
TSFC=(f./I_sp_a)*1000*3600;
TSFCmax = 24000*3600/((5000/1.0327)*65);
D_e = 2*sqrt((65000./I_sp_a)./(rhoe.*ve*pi));

for i = 1:length(Ttot2)
    if TSFC(i)>=TSFCmax
        TSFC(i) = NaN;
    end
end

[TSFC_min, F_min] = min(TSFC);

if ismember('plot', varargin)
    figure()
    plot(Ttot2, TSFC, 'LineWidth', 1);
    hold on;
    %yline(TSFC_min, '--k');
    xline(2100, 'Color', 'r', 'LineStyle','--')
    yline(TSFCmax, 'Color', 'r', 'LineStyle','-.')
    plot(700, TSFCmax, '>', 'LineWidth', 1, 'Color', 'r')
    ylim([180, 290])
    xlim([700, 2200])
    plot(2100, 180, '^', 'LineWidth', 1, 'Color', 'r')
    plot(Ttot2(F_min), TSFC_min, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor','b', 'MarkerSize', 5);
    lgd = legend("$TSFC$", "$T_{AB, max}$","$TSFC_{max}$", "", "", "$TSFC_{min}$", 'Interpreter','latex');
    lgd.FontSize = 8;
    lgd.IconColumnWidth = 18;
    xlabel("$T_{tot, 2}\,[K]$", "Interpreter","latex");
    ylabel("$TSFC \,\,\left[ \,\frac{kg}{kN\,h} \,\right] $", "Interpreter","latex")
    grid on

end

supCruise.m_a=T_supercruise/I_sp_a(F_min);
supCruise.f = f(F_min);
supCruise.I_sp_a=I_sp_a(F_min);
supCruise.TSFC=TSFC(F_min);
supCruise.M=M_supercruise;
supCruise.v0=v0_supercruise;

