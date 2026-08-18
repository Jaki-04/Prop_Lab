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

options = optimoptions('fsolve','Display','none');
options2 = optimoptions('fminunc','Display','none');

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

max_toPlot1 = zeros(1, length(cone_angles));
max_toPlot2 = zeros(1, length(cone_angles));
max_toPlot3 = zeros(1, length(cone_angles));

for id1 = 1:length(cone_angles)
        max_toPlot1(id1) = max(max(ptot_grid(id1, :, :)));
end

for id2 = 1:length(cone_angles)
        max_toPlot2(id2) = max(max(ptot_grid(:, id2, :)));
end

for id3 = 1:length(cone_angles)
        max_toPlot3(id3) = max(max(ptot_grid(:, :, id3)));
end


%PLOT DELLA PRESSIONE AL VARIARE DEGLI ANGOLI DI RAMPA
    % [max1, id1] = max(max_toPlot1);
    % [max2, id2] = max(max_toPlot2);
    % [max3, id3] = max(max_toPlot3);
    % figure()
    % tlrampa =tiledlayout(1,3,'TileSpacing','loose');
    % title(tlrampa,"" );
    % %plot spostamento x
    % nexttile
    % plot(rad2deg(cone_angles), max_toPlot1,'LineWidth',1, 'Color','blue', 'LineWidth',1)
    % hold on;
    % plot(rad2deg(cone_angles(id1)), max1, 'o', 'Color','k');
    % plot(rad2deg(cone_angles(1:id1)), max1*ones(size(ramp_angles1(1:id1))), 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
    % pbaspect([1,1,1]);
    % ylabel('$P_{tot}$ [Pa]', 'Interpreter','latex')
    % xlabel('Angolo della prima rampa [deg]', 'Interpreter','latex')
    % grid on;
    % text(rad2deg(cone_angles(2)), max2 - 0.000025, '0.73107');
    % ylim([0.7304, 0.7311])
    % 
    % nexttile
    % plot(rad2deg(ramp_angles1), max_toPlot2,'LineWidth',1, 'Color','#FF8800', 'LineWidth',1)
    % hold on;
    % plot(rad2deg(ramp_angles1(id2)), max2, 'o', 'Color','k');
    % plot(rad2deg(ramp_angles1(1:id2)), max2*ones(size(ramp_angles1(1:id2))), 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
    % pbaspect([1,1,1])
    % ylabel('$P_{tot}$ [Pa]', 'Interpreter','latex')
    % xlabel('Angolo della seconda rampa [deg]', 'Interpreter','latex')
    % text(rad2deg(ramp_angles1(2)), max2 - 0.000025, '0.73107');
    % grid on;
    % ylim([0.7304, 0.7311])
    % 
    % nexttile
    % plot(rad2deg(ramp_angles2), max_toPlot3,'LineWidth',1, 'Color','red', 'LineWidth',1)
    % hold on;
    % plot(rad2deg(ramp_angles2(id3)), max3, 'o', 'Color', 'k');
    % plot(rad2deg(ramp_angles2(1:id3)), max3*ones(size(ramp_angles1(1:id3))), 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
    % pbaspect([1,1,1])
    % ylabel('$P_{tot}$ [Pa]', 'Interpreter','latex')
    % xlabel('Angolo della terza rampa [deg]', 'Interpreter','latex')
    % text(rad2deg(ramp_angles2(2)), max2 - 0.000025, '0.73107');
    % grid on;
    % ylim([0.7304, 0.7311])

[ptot_ratio, linear_idx] = max(ptot_grid(:), [], 'omitnan');
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

%sprintf('Rendimento massimo della presa pi_d=%f con valori %f° %f° %f°', ptot_ratio, rad2deg(delta1), rad2deg(delta2), rad2deg(delta3))

% % ALTRO PLOT DEL RENDIMENTO
% % 1. Genera la griglia 3D corretta per gli assi cartesiani (X=colonne, Y=righe, Z=profondità)
% [X, Y, Z] = meshgrid(rad2deg(cone_angles), rad2deg(ramp_angles1), rad2deg(ramp_angles2));
% V = ptot_grid;
% 
% % 2. Definisci il valore centrale e la tolleranza del range
% %valore_centro = (0.731071+0.728733)/2; 
% %tolleranza = (0.731071-0.728733)/2; % Se il grafico è vuoto, aumenta questo valore (es. 0.1 o 0.2)
% 
% valore_min = 0.73;
% valore_max = 0.732;
% 
% % 3. Trova i punti nel range ed estrai le coordinate cartesiane
% maschera = (V >= valore_min) & (V <= valore_max);
% x_punti = X(maschera);
% y_punti = Y(maschera);
% z_punti = Z(maschera);
% valori_punti = V(maschera);
% 
% % 4. Controllo di sicurezza e Plot
% figure;
% if isempty(x_punti)
%     % Se non trova nulla, ti stampa a schermo i valori reali della tua matrice per capire dove sbatte
%     fprintf('ERRORE: Nessun punto trovato tra %g e %g.\n', valore_min, valore_max);
%     fprintf('I valori della tua ptot_grid vanno da un MIN di %g a un MAX di %g.\n', min(V(:)), max(V(:)));
%     error('Allarga la tolleranza o cambia il valore_centro.');
% else
%     % Disegna i punti. 40 è la dimensione del pallino, cambiala se serve
%     scatter3(x_punti, y_punti, z_punti, 20, valori_punti, 'filled', 'LineWidth',0.5);
% end
% 
% % 5. Impostazioni grafiche e nomi degli assi corretti
% view(3); 
% grid on; 
% axis tight;
% colorbar; % Mostra la barra con i valori esatti di ptot_grid relativi ai punti
% colormap('cool'); 
% hold on;
% plot3(11.8, 14.4, 17.7,'o', 'MarkerSize', 5, 'MarkerFaceColor','r')
% text(11.8, 14.4, 17.7,'\textbf{MAX}', 'Color', 'r', 'Interpreter','latex', 'FontSize',11, 'VerticalAlignment', 'cap')
% xlabel('$\delta_1$', 'Interpreter','latex'); 
% ylabel('$\delta_2$', 'Interpreter','latex'); 
% zlabel('$\delta_3$', 'Interpreter','latex');


%% Calcolo dimensioni della spina

% delta1 = 
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
    ptot1 = ptot_ratio*p*(1+M^2*(g-1)/2)^(g/(g-1));
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

% Nessun diffusore per il supersonico
T2=T1;
p2=p1;
rho2=rho1;
A2=Ain_sup;
M2=M_in;

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
phi_inf = deg2rad(7); 

G = Geometry();
rt_c = G.Rc_t;           % Dati dal compressore
rh_c = G.Rc_h1;          % Dati dal compressore

rin_sub=sqrt(rt_c^2-Ain_sub/pi);
deltar_inf=rin_sub-rh_c;
L_diff = deltar_inf/tan(phi_inf);
widthlength_ratio = L_diff/(rt_c-rin_sub);

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