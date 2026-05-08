function [pf, Tf, rhof, Mf, Af] = Presa(varargin)
%% Ottimizzazione della spina supersonica

Air = Air_parameters('25000');
p=Air.p;
T=Air.T;
g = Air.g; 
R = Air.R;

supCruise = Supercruise_parameters();

M = supCruise.M;        % Mach di ingresso alla presa
m_a_sup = supCruise.m_a;

% Angoli del cono e i due angoli di rampa su cui ciclare
cone_angles = deg2rad(11:0.1:12);
ramp_angles1 = deg2rad(14:0.1:15);
ramp_angles2 = deg2rad(17:0.1:18);

% Inizializzazione parametri
ptot_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
pstat_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
Mf_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
alpha1_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
alpha2_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
alpha3_grid = zeros(length(cone_angles), length(ramp_angles1), length(ramp_angles1));
alpha_max2=zeros(size(cone_angles));
alpha_max3=zeros(size(cone_angles));


options = optimoptions('fsolve','Display','none');
options2 = optimoptions('fminunc','Display','none');

% Funzioni utili

    % Funzione angolo di deflessione --> angolo dell'onda
    delta_fun = @(M,alpha) atan(2.*cot(alpha).*((M.^2.*sin(alpha).^2)-1)./(M.^2.*(g+cos(2.*alpha)) + 2));
    
    % Funzione Mach di uscita dall'onda
    M_exit = @(M, alpha, delta) sqrt( (1./sin(alpha-delta).^2) .* (1+M.^2.*sin(alpha).^2.*(g-1)/2) ./ (g.*M.^2.*sin(alpha).^2-(g-1)/2) );
    
    % Funzione perdita di pressione totale a cavallo dell'onda
    p_tot_ratio = @(Mn) ( ((g+1).*Mn.^2./(2+(g-1).*Mn.^2)).^(g./(g-1))).*((g+1)./(2.*g.*Mn.^2-g+1)).^(1/(g-1));

    % Funzione guadagno di pressione statica a cavallo dell'onda
    p_stat_ratio = @(Mn) 1+(2*g./(g+1)).*(Mn.^2-1);


% Ciclo sugli angoli del cono di ingresso


for ramp1_id = 1:length(ramp_angles1)
    
    ramp_angle1 = ramp_angles1(ramp1_id);

    for ramp2_id = 1:length(ramp_angles2)

        ramp_angle2 = ramp_angles2(ramp2_id);

        % Prima onda obliqua
        alpha1 = fsolve(@(alpha) delta_fun(M, alpha)-cone_angles, deg2rad(35)*ones(size(cone_angles)), options);
        M1=M_exit(M, alpha1, cone_angles);
        Mn_supercruise = M.*sin(alpha1);
        alpha_max1= fminunc(@(alpha)-delta_fun(M, alpha), deg2rad(70), options2);
        delta_max1 = delta_fun(M, alpha_max1);
        
        % Seconda onda obliqua
        alpha2 = fsolve(@(alpha) delta_fun(M1, alpha)-ramp_angle1, deg2rad(45)*ones(size(cone_angles)), options);
        M2=M_exit(M1, alpha2, ramp_angle1);
        Mn1 = M1.*sin(alpha2);
        for i = 1:length(cone_angles)
            alpha_max2(i)= fminunc(@(alpha)-delta_fun(M1(i), alpha), deg2rad(70), options2);
        end
        delta_max2 = delta_fun(M1, alpha_max2);
    
        % Terza onda obliqua
        alpha3 = fsolve(@(alpha) delta_fun(M2,alpha)-ramp_angle2, deg2rad(55)*ones(size(cone_angles)), options);
        M3=M_exit(M2, alpha3, ramp_angle2);
        Mn2 = M2.*sin(alpha3);
        alpha_max2 = zeros(size(M2));
        for i = 1:length(cone_angles)
           alpha_max3(i) = fminunc(@(alpha)-delta_fun(M2(i), alpha), deg2rad(90), options2);
        end
        delta_max3 = delta_fun(M2, alpha_max3);
    
        % Onda normale
        alpha4 = pi/2;
        M4 = M_exit(M3, alpha4, zeros(size(ramp_angle1)));
        Mn3 = M3.*sin(alpha4);
        
        % Perdita complessiva 
        p_tot_final_ratio = p_tot_ratio(Mn_supercruise).*p_tot_ratio(Mn1).*p_tot_ratio(Mn2).*p_tot_ratio(Mn3);
    
        p_stat_final_ratio = p_stat_ratio(Mn_supercruise).*p_stat_ratio(Mn1).*p_stat_ratio(Mn2).*p_stat_ratio(Mn3);
        
        % Troncamento risultati antifisici
        for i=1:length(p_tot_final_ratio)
            if Mn1(i)<1 || Mn2(i)<1 || Mn3(i)<1 || delta_max1<= cone_angles(i) || delta_max2(i) <= ramp_angle1 || delta_max3(i) <= ramp_angle2
                p_tot_final_ratio(i) = NaN;
            end
        end
        
        % Riempio la colonna di p_grid
        ptot_grid(:, ramp1_id, ramp2_id) = p_tot_final_ratio;
        pstat_grid(:, ramp1_id, ramp2_id) = p_stat_final_ratio;
        Mf_grid(:, ramp1_id, ramp2_id) = M4;
        alpha1_grid(:, ramp1_id, ramp2_id) = alpha1;
        alpha2_grid(:, ramp1_id, ramp2_id) = alpha2;
        alpha3_grid(:, ramp1_id, ramp2_id) = alpha3;
    end
end

[p_tot_ratio_max, linear_idx] = max(ptot_grid(:), [], 'omitnan');
[id1, id2, id3] = ind2sub(size(ptot_grid), linear_idx);

delta1= cone_angles(id1);
delta2 = ramp_angles1(id2)+delta1;
delta3 = ramp_angles2(id3)+delta2;

p_ratio = pstat_grid(id1, id2, id3);
M_in = Mf_grid(id1, id2, id3);
a1 = alpha1_grid(id1, id2, id3);
a2 = alpha2_grid(id1, id2, id3);
a3 = alpha3_grid(id1, id2, id3);

% Plot della curva del rendimento della presa
% figure()
%     s=surf(rad2deg(cone_angles), rad2deg(ramp_angle1), ptot_grid');
%     s.EdgeColor = 'none';
%     hold on
%     plot3(rad2deg(cone_fin), rad2deg(ramp_fin), p_tot_ratio_max, 'o', 'MarkerFaceColor', 'r')
%     xlabel('Cone Angle')
%     ylabel('Ramp angles')
%     zlabel('$\pi_d$', 'Interpreter','latex')

sprintf('Rendimento massimo della presa pi_d=%f con valori %f° %f° %f°', p_tot_ratio_max, rad2deg(delta1), rad2deg(delta2), rad2deg(delta3))

%% Calcolo dimensioni della spina

% d1 = deg2rad(11.8);
% d2 = deg2rad(14.4);
% d3 = deg2rad(17.7);
% M_in = 0.7191;
% a1 = 0.4546;
% a2 = 0.5804;
% a3 = 0.8035;
% p_ratio = 39.5125;

    % Quantità dopo l'onda normale
    T1 = T*(1+M^2*(g-1)/2)/(1+M_in^2*(g-1)/2);
    p1 = p*p_ratio;
    Ttot1 = T*(1+M^2*(g-1)/2);
    ptot1 = p_tot_ratio_max*p*(1+M^2*(g-1)/2)^(g/(g-1));
    rho1=p1/(R*T1);
    v1 = M_in*sqrt(g*R*T1);

Ain_sup = m_a_sup/ (rho1*v1);             % Area di ingresso del diffusore supersonico
ds = 0.01;                                % Distanza (parallela alla rampa 3) dell'urto normale dal bordo superiore

s3 = @(h) h/tan(a3);
h3 = @(h) s3(h)*sin(delta3);

s2 = @(h) ((s3(h)*tan(delta3-delta2) + h)*cos(delta3-delta2))*(1/tan(a2) - 1/tan(a3+delta3-delta2));
h2 = @(h) s2(h)*sin(delta2);

h_temp = @(h) (s3(h)*tan(delta3-delta2) + h)*cos(delta3-delta2);
s_temp = @(h) h_temp(h)/tan(a3+delta3-delta2);
s1 = @(h) ((s2(h)+s_temp(h))*tan(delta2-delta1) + h_temp(h))*cos(delta2-delta1)*(1/tan(a1) - 1/tan(a2+delta2-delta1));
h1 = @(h) s1(h)*sin(delta1);

r_shock = @(h) h1(h)+h2(h)+h3(h);
Rmin_fun = @(h) r_shock(h) + ds*sin(delta3);
Rmax_fun = @(h) r_shock(h) + h*cos(delta3) +ds*sin(delta3);
Area_fun = @(h) pi * h * (Rmax_fun(h) + Rmin_fun(h));

Hin_sup = fsolve(@(h)Area_fun(h)-Ain_sup, 0.4, options);
r_in = Rmin_fun(Hin_sup);
R_in = Rmax_fun(Hin_sup);
L_sup = (R_in-ds*sin(delta3))/tan(a1+delta1)+ds*cos(delta3);    % Distanza dalla punta al cowl lip

% Diffusore subsonico dopo l'onda normale
M2 = 0.25;                                  % Mach target all'ingresso del postbruciatore
eta_diff = 0.97;
pi_presa = ((1+eta_diff*M_in^2*(g-1)/2)/(1+M_in^2*(g-1)/2))^(g/(g-1));    % Perdita di pressione totale

T2 = Ttot1/(1+M2^2*(g-1)/2);
v2 = M2*sqrt(R*g*T2);
p2 = ptot1*pi_presa/((1+M2^2*(g-1)/2)^(g/(g-1)));
rho2 = p2/(R*T2);

Area_ratio_sup = rho1*v1/(rho2*v2);
A2 = Area_ratio_sup*Ain_sup; 

%% Presa subsonica

Air = Air_parameters('12000');
p_sub=Air.p;
T_sub=Air.T;
rho_sub=Air.rho;
g = Air.g; 
R = Air.R;

subCruise = Subcruise_parameters();
M_sub = subCruise.M;
v0_sub = subCruise.v0;
m_a_sub = subCruise.m_a;

Ttot1_sub = T_sub*(1+M_sub^2*(g-1)/2);
ptot1_sub = p_sub*(1+M_sub^2*(g-1)/2)^(g/(g-1));

% Sezione di ingresso subsonica
Ain_sub = m_a_sub/(rho_sub*v0_sub);
rmax_spina = sqrt(R_in^2-Ain_sub/pi);
h0_d_sub = R_in-rmax_spina;
phi_inf = deg2rad(10); 

rt_c = 0.506560138592840;           % Dati dal compressore
rh_c = 0.253280069296420;           % Dati dal compressore

deltar_inf = rmax_spina-rh_c;
deltar_sup = rt_c-R_in;
L_diff = deltar_inf/tan(phi_inf);
L_sub = L_diff+sin(acos(r_in/rmax_spina))*rmax_spina;
phi_sup = atan(deltar_sup/L_diff);
widthlength_ratio = L_diff/h0_d_sub;

% Diffusore subsonico
eta_diff_sub = 0.97;               
M2_sub = 0.5;               % Mach target all'ingresso del compressore

pi_presa_sub = ((1+eta_diff_sub*M_sub^2*(g-1)/2)/(1+M_sub^2*(g-1)/2))^(g/(g-1));
T2_sub = Ttot1_sub/(1+M2_sub^2*(g-1)/2);
v2_sub = M2_sub*sqrt(R*g*T2_sub);
p2_sub = ptot1_sub*pi_presa_sub/((1+M2_sub^2*(g-1)/2)^(g/(g-1)));

rho2_sub = p2_sub/(R*T2_sub);
Area_ratio_sub = rho_sub*v0_sub/(rho2_sub*v2_sub);
A2_sub = Area_ratio_sub*Ain_sub;               % Area di ingresso al compressore

if ismember('sup', varargin) || ismember('supersonic', varargin) || ismember('Sup', varargin) || ismember('Supersonic', varargin)
    pf = p2;
    Tf = T2;
    rhof = rho2;
    Mf = M2;
    Af = A2;
else
    pf = p2_sub;
    Tf = T2_sub;
    rhof = rho2_sub;
    Mf = M2_sub;
    Af = A2_sub;
end

end
