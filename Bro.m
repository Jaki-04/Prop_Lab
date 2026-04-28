clc
clearvars;

f=0:0.0001:1;

p_25000=0.054*101325; % Sbagliato
T_25000=216;
Msup = 3.5;
pi_AB=0.97;
cp_GC=1300;
gamma=1.4;
gamma_CG = 1.3;
R=287.01;
v0 = Msup*sqrt(gamma*R*T_25000);
DH=42000000;
eta=0.92;
T=65000;

pi_presa = 0.8;
Ttot1=T_25000*(1+(gamma-1)/2*Msup^2);
ptot1=p_25000*(1+(gamma-1)/2*Msup^2)^(gamma/(gamma-1))*pi_presa;

ptot2=pi_AB*ptot1;
Ttot2=(cp_GC*Ttot1+f*DH*eta)./((1+f)*cp_GC);

ve = sqrt( 2*cp_GC*Ttot2*( 1-(p_25000/ptot2).^(gamma_CG-1)/(gamma_CG)) );

m_a=T./((1+f).*ve-v0);

I_sp = T./m_a;
TSFC=m_a.*f./T;
m_f =f.*m_a;
plot(f, m_f, f, m_a, f, ve.*0.01);
legend("$\dot{m}_f$", "$\dot{m}_a$", "$ve$", 'Interpreter','latex');

min=m_f(1);
minf=f(1);
for i=1:max(size(f))
    if m_f(i)<min
        min = m_f(i);
        minf=f(i);
    end
end
min_vec=ones(10, 1)*minf;
hold on
%plot(min_vec, 0:9, '--r')

%%

%f=(0.002:0.0001:0.005)';
b=(2:0.1:20);

%fmat = f*ones(1, max(size(b)));
bmat = ones(max(size(f)), 1)*b;

fmat = 0.003;

m_f =2;
m_a=m_f/fmat;

cp_GC = 1155;
cp_a = 1000;
gamma=1.4;
gamma_GC = 1.33;
Msub = 0.85;
pi_presa = 0.8;
pi_camera = 0.97;
eta=0.97;
DH = 42000000;
etac =0.92;
etat=0.92;
T=12000;
R=287.01;

p_25000=0.054*101325; 
T_25000=0;
v0 = Msub*sqrt(gamma*R*T_25000);

Ttot1=T_25000*( 1+(gamma-1)/2*Msub^2 );
ptot1=p_25000*( 1+(gamma-1)/2*Msub^2 )^(gamma/(gamma-1))*pi_presa;

ptot2=ptot1.*(bmat);
Ttot2=(bmat).^((gamma-1)/gamma);

ptot3=pi_camera.*ptot2;
Ttot3=(cp_GC*Ttot2+fmat.*DH*eta).\( ((1+fmat).*cp_GC));


Ttot4 = Ttot3 - (1/(etac*etat)).*( cp_a./((1+fmat).*cp_GC) ).* (Ttot2-Ttot1);
ptot4 = 0.97*ptot3;

ve = sqrt( 2.*cp_GC.*Ttot2.*( 1-(p_25000./ptot2).^(gamma_GC-1)./(gamma_GC)) );

%m_a=T./((1+fmat).*ve-v0);
T =  ((1+fmat)*ve-v0)*m_a;
% I_sp = T./m_a;
% TSFC=m_a.*fmat./T;
% m_f =fmat.*m_a;
plot(bmat(1, :), T)
legend("$\dot{m}_f$", "$\dot{m}_a$", "$ve$", 'Interpreter','latex')