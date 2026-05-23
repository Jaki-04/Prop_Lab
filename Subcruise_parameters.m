%% Studio Parametrico - Regime subsonico


function subCruise = Subcruise_parameters(varargin)

% Rendimenti e parametri
pi_noAB=0.99;       % Fissato
pi_b=0.95;          % Fissato
pi_d = 0.98;        % Fissato
eta_presa = 0.98;     % Fissato (presa obiettivo)
eta_b=0.98;         % Fissato
eta_n = 0.98;
et=0.9;
ec=0.9;
H_f=43000000;       % Fissato
eta_m=0.92;                            % Massima SOT accettato

% Aria a 12km
Air = Air_parameters('12000');
p= Air.p; 
T= Air.T;
cp_a = Air.cp;
g_a = Air.g;
R_a = Air.R;

% Caratteristiche aria e GC
cp_GC=Air.cp_GC;
g_GC =Air.g_GC;

% Parametri: beta e temperatura di ingresso in turbina
Ttot3=(1100:10:1897);
b=(2:0.05:40);

Tt3mat = Ttot3'*ones(1, max(size(b)));
bmat = ones(max(size(Ttot3)), 1)*b;

% Regime subsonico
M_subsonic = 0.85;
T_subsonic=25000;
v0_subsonic = M_subsonic*sqrt(g_a*R_a*T);

% Presa
pi_presa = ((1+eta_presa*M_subsonic^2*(g_a-1)/2)/(1+M_subsonic^2*(g_a-1)/2))^(g_a/(g_a-1));
Ttot1=T*( 1+(g_a-1)/2*M_subsonic^2 );
ptot1=p*( 1+(g_a-1)/2*M_subsonic^2 )^(g_a/(g_a-1))*pi_presa;

% Compressore
ptot2=ptot1.*bmat;
Ttot2_id = Ttot1.*( bmat.^((g_a-1)/g_a) );
% Compressore reale:
eta_c = ( bmat.^( (g_a-1)./g_a ) - 1 )./( bmat.^( (g_a-1)./(g_a*ec) ) - 1 );
Ttot2 = Ttot1 + (Ttot2_id - Ttot1)./eta_c;

% Diffusore
ptot_diff = pi_d*ptot2;

% C.C.
ptot3=pi_b.*ptot_diff;
SOT = Tt3mat.*2./(g_GC+1);
eps_cool_NGV = ((10/550 * (SOT-1100))/100).*(SOT>1100) ;     % Spillamento per raffreddare la turbina (preso dal grafico)
eps_cool_blades = (12/350 * (SOT-1250)/100).*(SOT>1250);
eps_cool = eps_cool_NGV+eps_cool_blades;

fmat = ((1-eps_cool).*cp_a.*Ttot2-Tt3mat.*cp_GC+eps_cool.*Ttot2.*cp_GC)./(Tt3mat.*cp_GC-H_f.*eta_b);

% Turbina reale
Ttot4 = (Tt3mat.*cp_GC.*(1+fmat-eps_cool) + eps_cool.*cp_a.*Ttot2 - (cp_a./eta_m).* (Ttot2-Ttot1)) ./ ( (1+fmat-eps_cool).*cp_GC +eps_cool.*cp_a );
% Turbina ideale
tau_T= Ttot4./Tt3mat;
pi_t = tau_T.^( g_GC./( et.*(g_GC-1)) );
ptot4 = ptot3.*pi_t;

% Post bruciatore (spento)
ptot_AB = ptot4.*pi_noAB;
Ttot_AB = Ttot4;

% Ugello
T_ratio = (p./ptot_AB).^( (g_GC-1)./(g_GC) );
ve = sqrt( 2.*cp_GC.*Ttot_AB.*eta_n.*(1-T_ratio) );
I_sp_a = ((1+fmat).*ve-v0_subsonic);


% Parametri di merito
TSFC=fmat./I_sp_a;

for i=1:size(I_sp_a, 1)
    for j=1:size(I_sp_a, 2)
        if I_sp_a(i, j)<=0 || TSFC(i, j)<=0 || TSFC(i, j)>=1e-4 %|| eps_cool(i, j)>0.05
            I_sp_a(i, j)=NaN;
            TSFC(i, j) = NaN;
        end
    end
end

maxValue = max(max(I_sp_a));
[ind11, ind12] = find(I_sp_a==maxValue);
bImax=zeros(size(Ttot3));
I_sp_a_curve=zeros(size(Ttot3));
for temp=1:length(Ttot3)
    [Tmaxb, idmax] = max(I_sp_a(temp, :));
    bImax(temp) = b(idmax);
    I_sp_a_curve(temp) = I_sp_a(temp, idmax);
end
[target, Tind] = min(abs(I_sp_a_curve-0.95*maxValue));
bind=find(b==bImax(Tind));

if ismember('plot', varargin)
    figure()
    contourf(Ttot3, b, TSFC'*3600*1000, 10, "FaceAlpha",0.75);
    xlabel("$T_{tot, 4}\,[K]$", 'Interpreter','latex')
    c = colorbar;
    c.Label.String = "$TSFC\, \left[\, \frac{Kg}{kN\,h} \,\right]$";
    c.Label.Interpreter = 'latex';
     xline(1250*(g_GC+1)/2, 'k--', 'LineWidth', 1)
     xline(1100*(g_GC+1)/2, 'k--', 'LineWidth', 1)
    ylabel("$\beta_c$", 'Interpreter','latex')
    figure()
    
    contourf(Ttot3, b, I_sp_a',11,"FaceAlpha",0.75);
     xlabel("$T_{tot, 4}\,[K]$", 'Interpreter','latex')
    c = colorbar;
    c.Label.String = "$I_{sp,a}\,\left[\frac{m}{s}\right]$";
    c.Label.Interpreter = 'latex';
    ylabel("$\beta_c$", 'Interpreter','latex')
    hold on
    plot(Ttot3, bImax, '--r', 'LineWidth', 1.5)
    plot(Ttot3(Tind), b(bind), 'or', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    text(Ttot3(Tind), b(bind), "C", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter','latex', 'Color', 'r', 'FontSize', 12);
    xline(1250*(g_GC+1)/2, 'k--', 'LineWidth', 1)
    xline(1100*(g_GC+1)/2, 'k--', 'LineWidth', 1)
    plot(Ttot3(ind11), b(ind12), 'sr', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
    text(Ttot3(ind11), b(ind12), "M", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter','latex', 'Color', 'r', 'FontSize', 12);
    legend("","$max\,(\,I_{sp,a}\,)$", "Condizione di progetto","","", "Massimo Impulso", 'Interpreter','latex')
end

subCruise.m_a = T_subsonic/I_sp_a(Tind, bind);
subCruise.f=fmat(Tind, bind);
subCruise.b=bmat(Tind, bind);
subCruise.I_sp_a=I_sp_a(Tind, bind);
subCruise.TSFC=TSFC(Tind, bind);
subCruise.M=M_subsonic;
subCruise.v0=v0_subsonic;