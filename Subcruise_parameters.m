%% Studio Parametrico - Regime subsonico


function subCruise = Subcruise_parameters(varargin)

% Rendimenti e parametri
pi_noAB=0.95;       % Fissato
pi_b=0.98;          % Fissato
pi_d = 0.97;        % Fissato
eta_presa = 0.97;     % Fissato (presa obiettivo)
eta_b=0.98;         % Fissato
eta_n = 0.92;
et=0.9;
ec=0.9;
H_f=43000000;       % Fissato
eta_m=0.92;
SOT_max = 1200;                              % Massima SOT accettato
eps_cool = (9/500 * (SOT_max-1100))/100 * 2;     % Spillamento per raffreddare la turbina (preso dal grafico)

% Caratteristiche aria e GC
cp_GC=1243;
g_GC = 1.3;

% Aria a 12km
Air = Air_parameters('12000');
p= Air.p; 
T= Air.T;
cp_a = Air.cp;
g_a = Air.g;
R_a = Air.R;

f=(0.01:0.0001:0.04)';
b=(2:0.05:40);

fmat = f*ones(1, max(size(b)));
bmat = ones(max(size(f)), 1)*b;

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
Ttot3= ((1-eps_cool)*cp_a.*Ttot2 + fmat .* H_f.*eta_b)./((1+fmat-eps_cool).*cp_GC);

% Turbina reale
Ttot4 = (Ttot3.*cp_GC.*(1+fmat-eps_cool) + eps_cool*cp_a*Ttot2 - (cp_a/eta_m).* (Ttot2-Ttot1)) ./ ( (1+fmat-eps_cool).*cp_GC +eps_cool*cp_a );
% Turbina ideale
tau_T= Ttot4./Ttot3;
pi_t = tau_T.^( g_GC./( et.*(g_GC-1)) );
ptot4 = ptot3.*pi_t;

% Post bruciatore (spento)
ptot_AB = ptot4.*pi_noAB;
Ttot_AB = Ttot4;

% Ugello
T_ratio = 1-(p./ptot_AB).^( (g_GC-1)./(g_GC) );
ve = sqrt( 2.*cp_GC.*Ttot_AB.*eta_n.*T_ratio );
m_a=T_subsonic./((1+fmat).*ve-v0_subsonic);
m_f =fmat.*m_a;

% Parametri di merito
I_sp_a = T_subsonic./m_a;
TSFC=m_f./T_subsonic;

for i=1:size(I_sp_a, 1)
    for j=1:size(I_sp_a, 2)
        if I_sp_a(i, j)<=0 || Ttot3(i, j)*2/(g_GC+1)>SOT_max || TSFC(i, j)<=0 || TSFC(i, j)>=1e-4
            I_sp_a(i, j)=NaN;
            TSFC(i, j) = NaN;
        end
    end
end

maxValue = max(max(I_sp_a));
[ind12, ind22] = find(I_sp_a==maxValue);

if ismember('plot', varargin)
    figure()
    t=surf(f, b, I_sp_a');
    view([1, 1, 0.25])
    t.EdgeColor='none';
    pbaspect([1, 1, 1])
    hold on;
    plot3(f(ind12), b(ind22), I_sp_a(ind12, ind22),'o', 'Color', 'r', 'MarkerFaceColor','r')
    x=t.XData;
    y=t.YData;
    z=t.ZData;
    %%Create vectors out of surface's XData and YData
    x=x(:,1);
    y=y(1,:);
    %%Divide the lengths by the number of lines needed
    xnumlines = 50; % 10 lines
    ynumlines =50; % 10 partitions
    xspacing = round(length(x)/xnumlines);
    yspacing = round(length(y)/ynumlines);
    %%Plot the mesh lines 
    % Plotting lines in the X-Z plane
    hold on
    for i = 1:yspacing:length(y)
        Y1 = y(i)*ones(size(x)); % a constant vector
        Z1 = z(i,:);
        plot3(x,Y1,Z1,'-k');
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(x)
        X2 = x(i)*ones(size(y)); % a constant vector
        Z2 = z(:,i);
        plot3(X2,y,Z2,'-k');
    end

    figure()
    t=surf(f, b, TSFC');
    view([1, 1, 1])
    t.EdgeColor='none';
    pbaspect([1, 1, 1])
    hold on;
    plot3(f(ind12), b(ind22),TSFC(ind12, ind22),'o', 'Color', 'r', 'MarkerFaceColor','r')

    % TSFC diminuisce con beta; per ogni beta ha un ottimo su f
    % I_sp_a ha un ottimo a un certo beta e un certo f

    
end

subCruise.m_a = m_a(ind12, ind22);
subCruise.f=fmat(ind12, ind22);
subCruise.b=bmat(ind12, ind22);
subCruise.I_sp_a=I_sp_a(ind12, ind22);
subCruise.TSFC=TSFC(ind12, ind22);
subCruise.M=M_subsonic;
subCruise.v0=v0_subsonic;