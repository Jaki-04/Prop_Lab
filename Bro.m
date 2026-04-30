%% Presa d'aria supersonica
[m_a, f, I_sp_a, TSFC] = Supercruise_parameters('plot');

M_supercruise = 3.5;        % Mach di ingresso alla presa
gamma_a =1.4;               % gamma dell'aria

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
    delta_fun = @(M,alpha) atan(2.*cot(alpha).*((M.^2.*sin(alpha).^2)-1)./(M.^2.*(gamma_a+cos(2.*alpha)) + 2));
    
    % Funzione Mach di uscita dall'onda
    M_exit = @(M, alpha, delta) sqrt( (1./sin(alpha-delta).^2) .* (1+M.^2.*sin(alpha).^2.*(gamma_a-1)/2) ./ (gamma_a.*M.^2.*sin(alpha).^2-(gamma_a-1)/2) );
    
    % Funzione perdita di pressione totale a cavallo dell'onda
    p_tot_ratio = @(Mn) ( ((gamma_a+1).*Mn.^2./(2+(gamma_a-1).*Mn.^2)).^(gamma_a./(gamma_a-1))).*((gamma_a+1)./(2.*gamma_a.*Mn.^2-gamma_a+1)).^(1/(gamma_a-1));


% Ciclo sugli angoli del cono di ingresso

ramp_id_max=1;
cone_id_max=1;
cone_id = 1;
M1_fin = 0;
M2_fin = zeros(size(ramp_angles));
M3_fin = zeros(size(ramp_angles));
M4_fin = zeros(size(ramp_angles));
alpha1_fin = 0;
alpha2_fin = zeros(size(ramp_angles));
alpha3_fin = zeros(size(ramp_angles));

for delta_cone = cone_angles

    % Prima onda obliqua
    alpha1 = fsolve(@(alpha) delta_fun(M_supercruise, alpha)-delta_cone, deg2rad(35), options);
    M1=M_exit(M_supercruise, alpha1, delta_cone);
    Mn_supercruise = M_supercruise.*sin(alpha1);
    
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
H= sqrt( m_a(F_min)/ (pi*rho_25000*v0_supercruise)); 
L_sup = H/tan(cone_fin+alpha1_fin);
l3 = H/tan(cone_fin+ramp_fin+ramp_fin+alpha3_fin);
l2 = H/tan(cone_fin+ramp_fin+alpha2_fin) - l3;
l1 = L_sup-l3-l2;