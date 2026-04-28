clc
clearvars;

% Rendimenti e parametri
pi_AB=0.97;
pi_b=0.98;
pi_presa = 0.8;
eta_b=0.98;
eta_AB=0.92;
eta_n = 0.95;
H_f=42000000;
eta_t=0.92;
eta_c =0.92;


% Caratteristiche aria e GC
cp_GC=1150;
cp_a=1000;
gamma_a=1.4;
gamma_GC = 1.3;
R_a=287.01;

% Aria a 25km


%% Crociera supersonica
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
ve = sqrt( 2*cp_GC*Ttot2*T_ratio);

% Portata d'aria target e portata massica corrispondente
m_a=T_supercruise./((1+f).*ve-v0);
m_f =f.*m_a;

% Parametri di merito
I_sp = T_supercruise./m_a;
TSFC=m_f./T_supercruise;
figure()
plot(f, TSFC*3600);
legend("$TSFC$", "$I_sp}$", "$ve$", 'Interpreter','latex');


%% Regime subsonico
p_12000= 19267; % Sbagliato
T_12000= 216.65;

f=(0.0045:0.00001:0.01)';
b=(2:0.01:25);

fmat = f*ones(1, max(size(b)));
bmat = ones(max(size(f)), 1)*b;
%bmat=7;

% Regime subsonico
M_subsonic = 0.85;
T_subsonic=12000;
v0 = M_subsonic*sqrt(gamma_a*R_a*T_25000);

% Presa
Ttot1=T_25000*( 1+(gamma_a-1)/2*M_subsonic^2 );
ptot1=p_25000*( 1+(gamma_a-1)/2*M_subsonic^2 )^(gamma_a/(gamma_a-1))*pi_presa;

% Compressore
ptot2=ptot1.*(bmat);
Ttot2 = Ttot1.*(bmat).^((gamma_a-1)/gamma_a);
%compressore reale:
%T2_id = Ttot1*(Ptot2/Ptot1)^( (gamma_a - 1) / gamma_a)
%Ttot2 = Ttot1 + (T2_id - Ttot1)/eta_c

% C.C.
ptot3=pi_b.*ptot2;
Ttot3= (cp_a.*Ttot2 + f .* H_f.*eta_b)./((1+f).*cp_GC);

% Turbina
Ttot4 = Ttot3 - (1/(eta_c*eta_t)).*( cp_a./((1+fmat).*cp_GC) ).* (Ttot2-Ttot1);
%turbina reale:
%T4_id = Ttot3 - (Ttot4 - Ttot3)./eta_t
%ptot4 = ptot3.*(T4_id./Ttot3).^((gamma_a-1)/gamma_a);
ptot4 = ptot3.*(Ttot4./Ttot3).^(gamma_GC/(gamma_GC-1));

% Ugello
ve = sqrt( 2.*cp_GC.*Ttot2.*( 1-(p_25000./ptot4).^( (gamma_GC-1)./(gamma_GC) ) ) );
m_a=T_subsonic./((1+fmat).*ve-v0);
m_f =fmat.*m_a;

% Parametri di merito
I_sp = T_subsonic./m_a;
TSFC=m_f./T_subsonic;

[mins, min_id] = min(TSFC(:, 1:end));
f_min = f(min_id);

figure()
plot(b, f_min);

%plot(f, TSFC(:, 1))
figure()
%surf(TSFC)
plot(f, TSFC)
ylim([-0.001, 0.001])
%legend("$\dot{m}_f$", "$\dot{m}_a$", "$ve$", 'Interpreter','latex')