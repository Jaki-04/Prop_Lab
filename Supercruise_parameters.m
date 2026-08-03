%% Studio parametrico - Crociera supersonica

function supCruise = Supercruise_parameters(varargin)

% Rendimenti e parametri
pi_AB = 0.97;       % Fissato
pi_d_sup = 0.73;     % Fissato (presa obiettivo)
eta_AB=0.98;
eta_n = 1;
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
%D_e = 2*sqrt((65000./I_sp_a)./(rhoe.*ve*pi));

TSFC(I_sp_a <= 0) = NaN;
TSFC(f <= 0) = NaN;

for i = 1:length(Ttot2)
    if TSFC(i)>=TSFCmax
        TSFC(i) = NaN;
    end
end

[TSFC_min, T_min] = min(TSFC);

[mod, Tind] = min(abs(TSFC(T_min:end)-1.025*TSFC_min));
Tind=T_min+Tind;

if ismember('plot', varargin)
    figure()
    plot(Ttot2, TSFC, 'LineWidth', 1.5, 'Color',[0.8500, 0.3250, 0.0980]);
    hold on;
    %yline(TSFC_min, '--k');
    xline(2100, 'Color', 'r', 'LineStyle','--')
    yline(TSFCmax, 'Color', 'r', 'LineStyle','-.')
    plot(700, TSFCmax, '>', 'LineWidth', 1, 'Color', 'r')
    ylim([180, 290])
    xlim([700, 2200])
    plot(2100, 180, '^', 'LineWidth', 1, 'Color', 'r')
    plot(Ttot2(T_min), TSFC_min, 'o', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], 'MarkerSize', 5)
    text(Ttot2(T_min)-10, TSFC_min-3, 'm','FontSize',12,'FontWeight','bold','Color',[0.8500, 0.3250, 0.0980],'Interpreter','latex')
    xlabel("$T_{tot, 7}\,[K]$", "Interpreter","latex");
    ylabel("$TSFC \,\,\left[ \,\frac{kg}{kN\,h} \,\right] $", "Interpreter","latex")
    plot(Ttot2(Tind), TSFC(Tind), 'o', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], 'MarkerSize', 5)
    text(Ttot2(Tind)-10, TSFC(Tind)-4, 'C', 'FontSize',12,'FontWeight','bold','Color',[0.8500, 0.3250, 0.0980],'Interpreter','latex')
    grid on
    
    ax = gca;
ax.TickLabelInterpreter = 'latex';

% --- Asse X ---
xt = xticks;                    % tick di default
xt = unique([xt, 2100]);        % forza la presenza del tick a 2100
xticks(xt);
xtl = "$" + string(xt) + "$";                      % wrappa tutto in $...$
xtl(abs(xt - 2100) < 1e-6) = "$T_{AB,\mathrm{max}}$";
xticklabels(xtl);

% --- Asse Y ---
yt = yticks;
yt = unique([yt, TSFCmax]);     % forza la presenza del tick a TSFCmax
yticks(yt);
ytl = "$" + string(yt) + "$";
ytl(abs(yt - TSFCmax) < 1e-6) = "$TSFC_{\mathrm{max}}$";
yticklabels(ytl);

figure()
plot(Ttot2, I_sp_a,'Linewidth', 1.5)
xline(2100, 'Color', 'r', 'LineStyle','--')

ax = gca;
ax.TickLabelInterpreter = 'latex';
% --- Asse X ---
xt = xticks;                    % tick di default
xt = unique([xt, 2100]);        % forza la presenza del tick a 2100
xticks(xt);
xtl = "$" + string(xt) + "$";                      % wrappa tutto in $...$
xtl(abs(xt - 2100) < 1e-6) = "$T_{AB,\mathrm{max}}$";
xticklabels(xtl);
grid on
hold on
plot(2100, 0, '^', 'LineWidth', 1, 'Color', 'r')
ylim([0,900])
xlabel("$T_{tot, 7}\,[K]$", "Interpreter","latex");
   ylabel("$I_{sp,a} \,\,\left[ \,\frac{m}{s} \,\right] $", "Interpreter","latex")
end

supCruise.m_a=T_supercruise/I_sp_a(Tind);
supCruise.f = f(Tind);
supCruise.I_sp_a=I_sp_a(Tind);
supCruise.TSFC=TSFC(Tind);
supCruise.M=M_supercruise;
supCruise.v0=v0_supercruise;

