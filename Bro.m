clc
clearvars;

% Rendimenti e parametri
pi_noAB=0.95;       % Fissato
pi_AB = 0.98;       % Fissato
pi_b=0.98;          % Fissato
pi_d = 0.97;        % Fissato
pi_presa = 0.8;     % Fissato (presa obiettivo)
eta_b=0.98;         % Fissato
eta_AB=0.92;
eta_n = 0.92;
et=0.9;
ec=0.9;
H_f=43000000;       % Fissato
eta_m=0.92;


% Caratteristiche aria e GC
cp_GC=1150;
cp_a=1000;
gamma_a=1.4;
gamma_GC = 1.33;
R_a=287.01;


%% Crociera supersonica

% Aria a 25km
p_25000= 2549; 
T_25000= 216.65;

f=0.0045:0.0001:0.1;

% Crociera supersonica
M_supercruise = 3.5;
v0 = M_supercruise*sqrt(gamma_a*R_a*T_25000);
T_supercruise=65000;

% Presa
Ttot1=T_25000*(1+(gamma_a-1)/2*M_supercruise^2);
ptot1=p_25000*(1+(gamma_a-1)/2*M_supercruise^2)^(gamma_a/(gamma_a-1))*pi_presa;

% AB
ptot2=pi_AB*ptot1;
Ttot2=(cp_a*Ttot1+f*H_f*eta_AB)./((1+f)*cp_GC);

% Ugello
T_ratio = 1-(p_25000/ptot2).^( (gamma_GC-1)/(gamma_GC) );
ve = sqrt( 2.*cp_GC.*Ttot2.*eta_n.*T_ratio);

% Portata d'aria target e portata massica corrispondente
m_a=T_supercruise./((1+f).*ve-v0);
m_f =f.*m_a;

% Parametri di merito

TSFC=m_f./T_supercruise;
figure()
plot(f, TSFC*3600);
ylim([0.001, 0.5])
legend("$TSFC$", "$I_sp}$", "$ve$", 'Interpreter','latex');


%% Regime subsonico

% Aria a 12km
p_12000= 19267; 
T_12000= 216.65;

f=(0.01:0.0001:0.04)';
b=(2:0.05:40);

fmat = f*ones(1, max(size(b)));
bmat = ones(max(size(f)), 1)*b;

% Regime subsonico
M_subsonic = 0.85;
T_subsonic=12000;
v0 = M_subsonic*sqrt(gamma_a*R_a*T_12000);

% Presa
Ttot1=T_12000*( 1+(gamma_a-1)/2*M_subsonic^2 );
ptot1=p_12000*( 1+(gamma_a-1)/2*M_subsonic^2 )^(gamma_a/(gamma_a-1))*pi_presa;

% Compressore
ptot2=ptot1.*bmat;
Ttot2_id = Ttot1.*( bmat.^((gamma_a-1)/gamma_a) );
% Compressore reale:
eta_c = ( bmat.^( (gamma_a-1)./gamma_a ) - 1 )./( bmat.^( (gamma_a-1)./(gamma_a*ec) ) - 1 );
Ttot2 = Ttot1 + (Ttot2_id - Ttot1)./eta_c;

% Diffusore
ptot_diff = pi_d*ptot2;

% C.C.
ptot3=pi_b.*ptot_diff;
Ttot3= (cp_a.*Ttot2 + f .* H_f.*eta_b)./((1+f).*cp_GC);

% Turbina reale
Ttot4 = Ttot3 - (1/eta_m).*( cp_a./((1+fmat).*cp_GC) ).* (Ttot2-Ttot1);
% Turbina ideale
tau_T= Ttot4./Ttot3;
pi_t = tau_T.^( gamma_GC./( et.*(gamma_GC-1)) );
ptot4 = ptot3.*pi_t;

% Post bruciatore (spento)
ptot_AB = ptot4.*pi_noAB;
Ttot_AB = Ttot4;

% Ugello
T_ratio = 1-(p_12000./ptot_AB).^( (gamma_GC-1)./(gamma_GC) );
ve = sqrt( 2.*cp_GC.*Ttot_AB.*eta_n.*T_ratio );
m_a=T_subsonic./((1+fmat).*ve-v0);
m_f =fmat.*m_a;

% Parametri di merito
I_sp_a = T_subsonic./m_a;
TSFC=m_f./T_subsonic;

% Ricerca del minimo

% f_min = zeros(1, length(b));
% b_fmin=zeros(1, length(b));
% fmin_id = zeros(1, length(b));
% bmin_id=zeros(1, length(b));

% for i = 1:max(size(b))
%     [maximum, max_id] = max(TSFC(1:end, i));
%     if max_id>=max(size(f))
%         break
%     else
%         [minimum, min_id] = min(abs(TSFC(max_id:end, i)));
%         f_min(i) = f(min_id+max_id-1);
%         fmin_id(i) = min_id;
%         bmin_id(i) = i;
%         b_fmin(i) = b(i);
%     end
% end
% 
% figure()
% plot(b_fmin, f_min);

% figure()
% hold on
% t=surf(f, b, I_sp');
% t.EdgeColor='none';

% TSFCmin=zeros(1, length(b));
% for i=1:length(fmin_id)
%     TSFCmin(i) = TSFC(fmin_id(i), bmin_id(i));
% end
% plot3(f_min, b_fmin, TSFCmin', '--k')
% pbaspect([1 2 3])
% %plot(f, TSFC)
% zlim([0, 0.0001])

for i=1:size(TSFC, 1)
    for j=1:size(TSFC, 2)
        if TSFC(i, j)<=0 || TSFC(i, j)>=1e-4 || Ttot3(i, j)>1400
            TSFC(i, j)=NaN;
        end
    end
end

minValue = min(min(TSFC));
[ind11, ind21] = find(TSFC==minValue);

for i=1:size(I_sp_a, 1)
    for j=1:size(I_sp_a, 2)
        if I_sp_a(i, j)<=0 || Ttot3(i, j)>1400
            I_sp_a(i, j)=NaN;
        end
    end
end

maxValue = max(max(I_sp_a));
[ind12, ind22] = find(I_sp_a==maxValue);

for i=1:size(I_sp_f, 1)
    for j=1:size(I_sp_f, 2)
        if I_sp_f(i, j)<=0 || Ttot3(i, j)>1400
            I_sp_f(i, j)=NaN;
        end
    end
end

maxValue = max(max(I_sp_f));
[ind13, ind23] = find(I_sp_f==maxValue);

figure()
    t=surf(f, b, I_sp_a');
    t.EdgeColor='none';
    pbaspect([1, 1, 1])
    hold on;
    plot(f(ind12), b(ind22),'o')
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
    t.EdgeColor='none';
    pbaspect([1, 1, 2])
    hold on;
    plot3(f(ind12), b(ind22),TSFC(ind12, ind22),'o')

% TSFC diminuisce con beta; per ogni beta ha un ottimo su f
% I_sp_a ha un ottimo a un certo beta e un certo f
