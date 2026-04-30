%% Presa d'aria supersonica

Air = Air_parameters('25000');
p=Air.p;
T=Air.T;
rho=Air.rho;
g = Air.g; 
R = Air.R;

supCruise = Supercruise_parameters();

M = supCruise.M;        % Mach di ingresso alla presa
v0 = supCruise.v0;

% Angoli del cono e i due angoli di rampa (uguali) su cui ciclare
cone_angles = deg2rad(5:0.1:25);
ramp_angles = deg2rad(5:0.1:25);

% Inizializzazione parametri
p_grid = zeros(length(cone_angles), length(ramp_angles));
p_tot_ratio_max =0;
idx=1;
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

ramp_id_max=1;
cone_id_max=1;
cone_id = 1;

for delta_cone = cone_angles

    % Prima onda obliqua
    alpha1 = fsolve(@(alpha) delta_fun(M, alpha)-delta_cone, deg2rad(35), options);
    M1=M_exit(M, alpha1, delta_cone);
    Mn_supercruise = M.*sin(alpha1);
    
    % Seconda onda obliqua
    alpha2 = fsolve(@(alpha) delta_fun(M1, alpha)-ramp_angles, deg2rad(45)*ones(size(ramp_angles)), options);
    M2=M_exit(M1, alpha2, ramp_angles);
    Mn1 = M1.*sin(alpha2);
   alpha_max1= fminunc(@(alpha)-delta_fun(M1, alpha), deg2rad(70), options2);
   delta_max1 = delta_fun(M1, alpha_max1);


    % Terza onda obliqua
    alpha3 = fsolve(@(alpha) delta_fun(M2,alpha)-ramp_angles, deg2rad(55)*ones(size(ramp_angles)), options);
    M3=M_exit(M2, alpha3, ramp_angles);
    Mn2 = M2.*sin(alpha3);
    alpha_max2 = zeros(size(M2));
    for i = 1:length(M2)
       [alpha_max2(i)] = fminunc(@(alpha)-delta_fun(M2(i), alpha), deg2rad(90), options2);
    end
    delta_max2 = delta_fun(M2, alpha_max2);

    % Onda normale
    alpha4 = pi/2;
    M4 = M_exit(M3, alpha4, zeros(size(ramp_angles)));
    Mn3 = M3.*sin(alpha4);
    
    % Perdita complessiva 
    p_tot_final_ratio = p_tot_ratio(Mn_supercruise).*p_tot_ratio(Mn1).*p_tot_ratio(Mn2).*p_tot_ratio(Mn3);

    p_stat_final_ratio = p_stat_ratio(Mn_supercruise).*p_stat_ratio(Mn1).*p_stat_ratio(Mn2).*p_stat_ratio(Mn3);
    
    % Troncamento risultati antifisici
    for i=1:length(p_tot_final_ratio)
        if Mn1(i)<1 || Mn2(i)<1 || Mn3(i)<1 || delta_max1<= ramp_angles(i) || delta_max2(i) <= ramp_angles(i)
            p_tot_final_ratio(i) = NaN;
        end
    end
    
    % Riempio la colonna di p_grid
    p_grid(idx, :) = p_tot_final_ratio;
    idx=idx+1;
    
    % Ricerca ottimo dell'iterazione corrente
    [currmax, ramp_id] = max(p_tot_final_ratio);
    if currmax > p_tot_ratio_max
        p_tot_ratio_max = currmax;
        ramp_fin = ramp_angles(ramp_id);
        cone_fin = delta_cone;
        ramp_id_max = ramp_id;
        cone_id_max = cone_id;
        M1_fin = M1;
        M2_fin = M2(ramp_id);
        M3_fin = M3(ramp_id);
        M4_fin = M4(ramp_id);
        alpha1_fin = alpha1;
        alpha2_fin = alpha2(ramp_id);
        alpha3_fin = alpha3(ramp_id);
        p_ratio_fin = p_stat_final_ratio(ramp_id);
    end

    if delta_cone== deg2rad(18)
        MS=M2;
    end

    cone_id=cone_id+1;
end

% Plot della curva del rendimento della presa
figure()
    s=surf(rad2deg(cone_angles), rad2deg(ramp_angles), p_grid');
    s.EdgeColor = 'none';
    hold on
    plot3(rad2deg(cone_fin), rad2deg(ramp_fin), p_tot_ratio_max, 'o', 'MarkerFaceColor', 'r')
    xlabel('Cone Angle')
    ylabel('Ramp angles')
    zlabel('$\pi_d$', 'Interpreter','latex')

sprintf('Rendimento massimo della presa pi_d=%f con valori delta_cone=%f°, delta_ramp=%f° ', p_tot_ratio_max, rad2deg(cone_fin), rad2deg(ramp_fin))

%% Calcolo dimensioni della presa

H= sqrt( supCruise.m_a/ (pi*rho*v0));                   % Raggio della presa
L_sup = H/tan(cone_fin+alpha1_fin);
l3 = H/tan(cone_fin+ramp_fin+ramp_fin+alpha3_fin);
l2 = H/tan(cone_fin+ramp_fin+alpha2_fin) - l3;
l1 = L_sup-l3-l2;

T_f = T*(1+M^2*(g-1)/2)/(1+M4_fin^2*(g-1)/2);
p_f = p*p_ratio_fin;
rho_f=p_f/(R*T_f);
v0_f = M4_fin*sqrt(g*R*T_f);

delta_anulo = cone_fin+2*ramp_fin;
A_anulo= supCruise.m_a/ (rho_f*v0_f);         % Area dell'anulo
fun = @(l) pi*(H^2*l-H/tan(delta_anulo)*l^2+1/(3*tan(delta_anulo)^2)*l^3)-A_anulo;
zf = fsolve(fun, 0.4);
rf = H-zf/(tan(delta_anulo));
h_anulo=(H-rf)/cos(delta_anulo);        % Dovrebbe venire 0.42286
h_anulo = 0.42286;

%% Presa subsonica

h_anulo = 0.42286;
Air = Air_parameters('12000');
p=Air.p;
T=Air.T;
rho=Air.rho;
g = Air.g; 
R = Air.R;

subCruise = Subcruise_parameters();
M = subCruise.M;
v0 = subCruise.v0;
m_a = subCruise.m_a;

% Sezione di ingresso subsonica
A_anulo_sub = m_a/(rho*v0);
R_top = sqrt(H^2-A_anulo_sub/pi);     % Raggio necessario alla sommità della spina per avere l'area richiesta
h0_duct = H-R_top;     % Altezza iniziale della strozzatura di ingresso in subsonico
fun = @(l) pi*(H^2*l-H/tan(delta_anulo)*l^2+1/(3*tan(delta_anulo)^2)*l^3)-A_anulo_sub;
zf = fsolve(fun, 0.4);
rf = H-zf/(tan(delta_anulo));
h_anulo=(H-rf)/cos(delta_anulo); 