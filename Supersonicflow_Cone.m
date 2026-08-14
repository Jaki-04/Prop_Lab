clear all
clc

k = 1.4; 
M_inf = 3.5;
format long
% alpha1_vec = deg2rad(25:0.2:27);
% alpha2_vec = deg2rad(32:0.2:35);
% alpha3_vec = deg2rad(43:0.2:45);
alpha1_vec = deg2rad(25.4);
alpha2_vec = deg2rad(32.4);
alpha3_vec = deg2rad(44);
p_ratio_low = 1;
p_ratio_high = 1;
ptot_grid_low = ones(length(alpha1_vec), length(alpha2_vec), length(alpha3_vec));
ptot_grid_high = ones(length(alpha1_vec), length(alpha2_vec), length(alpha3_vec));

for a1_id = 1:length(alpha1_vec)

    alpha1 = alpha1_vec(a1_id);

    for a2_id = 1:length(alpha2_vec)
    
        alpha2 = alpha2_vec(a2_id);

        for a3_id = 1:length(alpha3_vec)

            alpha3 = alpha3_vec(a3_id);
            alpha_it = [alpha1, alpha2, alpha3, deg2rad(90)];
            M_curr_low = M_inf;
            M_curr_high = M_inf;
            change_angle = 0;

            for i=1:4

                alpha_curr = alpha_it(i);
                
                if i<4
                    % Cone Angle Calculated from Theta-Beta-M function-file
                    [theta_c_r] = Theta_Beta_M_V2(0,alpha_curr,M_curr_low,k,'rad');                      
                    % Using oblique shock relations to get M behind shock
                    Mn_curr = M_curr_low*sin(alpha_curr);                                                % Normal Mach number pre-shock, Eqn. 4.7 (pg. 135)
                    Mn2 = sqrt((1 + ((k-1)/2)*Mn_curr^2)/((k*Mn_curr^2) - ((k-1)/2)));              % Normal Mach number post-shock, Eqn. 3.51 (pg. 89)
                    M2  = Mn2/(sin(alpha_curr-theta_c_r));                                   % Post-shock Mach number, Eqn. 4.12 (pg. 135)
                else
                    P0_P01 = ((((k+1)*M_curr_low^2)/((k-1)*M_curr_low^2+2))^(k/(k-1))) * ...
                    (((k+1)/(2*k*M_curr_low^2-(k-1)))^(1/(k-1)));
                    Mf_low = sqrt((1+M_curr_low.^2.*(k-1)/2) ./ (k.*M_curr_low.^2-(k-1)/2) );
                    p_ratio_low = p_ratio_low*(1+(2*k./(k+1)).*(M_curr_low.^2-1));
                    break
                end

                if M2>=M_curr_low || Mn2>=1 || (i<4 && M2<1)
                    ptot_grid_low(a1_id, a2_id, a3_id) = NaN;
                    change_angle = 1;

                    break;
                end
    
                    % Initial conditions
                V      = ((2/((k-1)*M2^2))+1)^(-1/2);                                   % Total velocity    [m/s]
                Vtheta = -V*(sin(alpha_curr-theta_c_r));                                 % Angular velocity  [m/s]
                Vr     = V*(cos(alpha_curr-theta_c_r));                                  % Radial velocity   [m/s]
                
                    %Theta range and initial conditions of Vr and Vr'
                thetarange = [alpha_curr; 1e-3];                                       % Integrate from shock angle to ~ 0 degrees
                V_init    = [Vr; Vtheta];                                               % Initial values for Vr and dVr/dTheta (Vtheta)
               
                options = odeset('Events',@EVENTS, 'Reltol', 2.22045e-10);
                sol     = ode23t(@TM_Equations,thetarange,V_init,options,k);            
                
                    % If the method worked, calculate angle and Mach number at the cone
                [n,m]  = size(sol.y);                                                   % Get the size of the solution array
                if (n > 0 && m > 0 && isempty(sol.ie) ~= 1)                             % If the solution converged
                    theta_c_r = sol.xe;                                                 % Final cone angle [rad]
                    theta_c_d = theta_c_r.*(180/pi);                                    % Final cone angle [deg]
                    Vc2       = sol.ye(1)^2 + sol.ye(2)^2;                              % Total velocity squared [m/s]^2
                    Mc        = ((1.0/Vc2-1)*(k-1)/2)^-0.5;                             % Mach number at cone surface
            
                else
                    ptot_grid_low(a1_id, a2_id, a3_id) = NaN;
                    change_angle = 1;

                    break
                end
                
                
                    % Oblique Shock Relations
                rho2_rho1 = (k+1)*Mn_curr^2/((k-1)*Mn_curr^2+2);                                % Eqn. 4.8 (pg. 135)
                P2_P1     = 1 + 2*k*(Mn_curr^2-1)/(k+1);                                    % Eqn. 4.9 (pg. 135)
                P02_P     = (1+(k-1)*Mc.^2/2).^(k/(k-1));                           % Eqn. 3.30 (pg. 80)
                P02_P2     = (1+(k-1)*M2^2/2)^(k/(k-1));                                % Eqn. 3.30 (pg. 80)
            
                P_P1 = P02_P2*P2_P1./P02_P;
                P0_P01 = P_P1*((1+(k-1)*Mc^2/2)/(1+(k-1)*M_curr_low^2/2))^(k/(k-1));
            
                ptot_grid_low(a1_id, a2_id, a3_id) = ptot_grid_low(a1_id, a2_id, a3_id)*P0_P01;
                M_curr_low = M2;
            end
            
            for i=1:4

                alpha_curr = alpha_it(i);

                if i<4
                    % Cone Angle Calculated from Theta-Beta-M function-file
                    [theta_c_r] = Theta_Beta_M_V2(0,alpha_curr,M_curr_high,k,'rad');                      
                    %f = rad2deg(theta_c_r)
                    % Using oblique shock relations to get M behind shock
                    Mn_curr = M_curr_high*sin(alpha_curr);                                                % Normal Mach number pre-shock, Eqn. 4.7 (pg. 135)
                    Mn2 = sqrt((1 + ((k-1)/2)*Mn_curr^2)/((k*Mn_curr^2) - ((k-1)/2)));              % Normal Mach number post-shock, Eqn. 3.51 (pg. 89)
                    M2  = Mn2/(sin(alpha_curr-theta_c_r));                                   % Post-shock Mach number, Eqn. 4.12 (pg. 135)
                else
                    P0_P01 = ((((k+1)*M_curr_high^2)/((k-1)*M_curr_high^2+2))^(k/(k-1))) * ...
                    (((k+1)/(2*k*M_curr_high^2-(k-1)))^(1/(k-1)));
                    Mf_high = sqrt((1+M_curr_high.^2.*(k-1)/2) ./ (k.*M_curr_high.^2-(k-1)/2) );
                    ptot_grid_high(a1_id, a2_id, a3_id) = ptot_grid_high(a1_id, a2_id, a3_id) * P0_P01;
                    break
                end

                if M2>=M_curr_high || Mn2>=1 || (i<4 && M2<1)
                    ptot_grid_high(a1_id, a2_id, a3_id) = NaN;
                    change_angle = 1;

                    break;
                end
    
                    % Initial conditions
                V      = ((2/((k-1)*M2^2))+1)^(-1/2);                                   % Total velocity    [m/s]
                Vtheta = -V*(sin(alpha_curr-theta_c_r));                                 % Angular velocity  [m/s]
                Vr     = V*(cos(alpha_curr-theta_c_r));                                  % Radial velocity   [m/s]
                
                    %Theta range and initial conditions of Vr and Vr'
                thetarange = [alpha_curr; 1e-3];                                       % Integrate from shock angle to ~ 0 degrees
                V_init    = [Vr; Vtheta];                                               % Initial values for Vr and dVr/dTheta (Vtheta)
               
                options = odeset('Events',@EVENTS, 'Reltol', 2.22045e-10);
                sol     = ode23t(@TM_Equations,thetarange,V_init,options,k);            
                
                    % If the method worked, calculate angle and Mach number at the cone
                [n,m]  = size(sol.y);                                                   % Get the size of the solution array
                if (n > 0 && m > 0 && isempty(sol.ie) ~= 1)                             % If the solution converged
                    theta_c_r = sol.xe;                                                 % Final cone angle [rad]
                    theta_c_d = theta_c_r.*(180/pi);                                    % Final cone angle [deg]
                    Vc2       = sol.ye(1)^2 + sol.ye(2)^2;                              % Total velocity squared [m/s]^2
                    Mc        = ((1.0/Vc2-1)*(k-1)/2)^-0.5;                             % Mach number at cone surface
            
                else
                    ptot_grid_high(a1_id, a2_id, a3_id) = NaN;
                    change_angle = 1;

                    break
                end
                
                
                    % Oblique Shock Relations
                rho2_rho1 = (k+1)*Mn_curr^2/((k-1)*Mn_curr^2+2);                                % Eqn. 4.8 (pg. 135)
                P2_P1     = 1 + 2*k*(Mn_curr^2-1)/(k+1);                                    % Eqn. 4.9 (pg. 135)
                P02_P     = (1+(k-1)*Mc.^2/2).^(k/(k-1));                           % Eqn. 3.30 (pg. 80)
                P02_P2     = (1+(k-1)*M2^2/2)^(k/(k-1));                                % Eqn. 3.30 (pg. 80)
            
                P_P1 = P02_P2*P2_P1./P02_P;
                P0_P01 = P_P1*((1+(k-1)*Mc^2/2)/(1+(k-1)*M_curr_high^2/2))^(k/(k-1));
            
                ptot_grid_high(a1_id, a2_id, a3_id) = ptot_grid_high(a1_id, a2_id, a3_id)*P0_P01;
                M_curr_high = Mc;
            end
            

            if change_angle
                continue
            end

        end
    end
end

function [val,last,dir] = EVENTS(theta,z,gam)

% Event Handling for ode15s Solver
% PURPOSE
% - Stops the ODE solution when the cone wall is reached
% INPUTS
% - theta : Integration angle [rad]
% - z     : Angular and radial velocities
% - gam   : Ratio of specific heats
% 
% OUTPUTS
% - val      : Value of the ith event function
%              Event is triggered when value is 0
% - last     : Integration terminates at a zero of the event function
%                  if this value is 1, otherwise it's 0
% - dir      : If all zeros are to be located, 0
%                If zeros where event fnc decreasing, -1
%                If zeros where event fnc increasing, +1

val  = z(2);   % evento triggerato quando Vθ attraversa zero
last = 1;      % ferma all'primo zero
dir  = 1;      % Vθ parte negativa e sale verso zero (direzione crescente)

% Trigger when angular velocity switches sign
% - Angular velocity needs to be zero at the cone surface (solid surface)

if (z(2) > 0)
    val = 0;              % If angular velocity z(2) crosses from
end                       %  being negative to positive, trigger event

end


function [sol] = Theta_Beta_M_V2(thetaInf,betaInf,MInf,k,degrad)

% PURPOSE
% Use the Theta-Beta-M relation to obtain the flow deflection angle from
%  the shock wave angle and freestream Mach number

% REFERENCES
% Modern Compressible Flow, Anderson, pg. 136, eqn. 4.17

% INPUTS
% thetaInf : Cone angle to horizontal [deg/rad] 
% betaInf  : Shockwave angle to horizontal [deg/rad]
% MInf     : Freestream Mach number 
% k        : Ratio of specific heats 
% degrad   : Specify whether calculations are in degrees or radians [str]

% OUTPUTS
% sol      : Solution value, depending on what is being solved for

% Convert between degrees and radians if necessary
if (strcmp(degrad,'deg'))
    thetaInf = thetaInf*(pi/180);
    betaInf  = betaInf*(pi/180);
elseif (strcmp(degrad,'rad'))
    % Do nothing, computations are in radians
end

% -------------------- SOLVING FOR: THETA ---------------------------------
if (thetaInf == 0)
    
    % Check against Mach wave angle equation
    if (betaInf <= asin(1/MInf))
        sol = 0;
        return;
    end
    
    % Numerator and denominator
    N = ((MInf*sin(betaInf))^2)-1;
    D = (MInf^2*(k+cos(2*betaInf)))+2;
    
    % Wedge angle from the Theta-Beta-M relation
    theta = atan(2*cot(betaInf)*(N/D));
    
    % Set solution variable based on input deg/rad
    if (strcmp(degrad,'deg'))                                               % If the output is in degrees
        sol = theta*(180/pi);
    elseif (strcmp(degrad,'rad'))                                           % If the output is in radians
        sol = theta;
    end
end
end
%{ 
----------------- SOLVING FOR: MACH NUMBER ------------------------------
elseif (MInf == 0)
    
    % Find zero of theta-beta-M equation
    myfun = @(t,b,g,M) (2*cot(b)*((((M*sin(b))^2)-1)/...
                        ((M^2*(g+cos(2*b)))+2))) - tan(t);
    t   = thetaInf;                                                         % Given cone angle [rad]
    b   = betaInf;                                                          % Given shock angle [rad]
    g   = k;                                                                % Given ratio of specific heats []
    fun = @(M) myfun(t,b,g,M);                                              % Set given variables in myfun
    M   = fzero(fun,1);                                                     % Solve for M, starting at M = 1
    
    % Set solution variable
    sol = M;                                                                % Set the solution variable
    
% -------------------- SOLVING FOR: BETA ----------------------------------
elseif (betaInf == 0)
    
    % Find zero of theta-beta-M equation
    myfun = @(t,b,g,M) (2*cot(b)*((((M*sin(b))^2)-1)/...
                        ((M^2*(g+cos(2*b)))+2))) - tan(t);
    t   = thetaInf;                                                         % Given cone angle [rad]
    M   = MInf;                                                             % Given Mach number [rad]
    g   = k;                                                              % Given ratio of specific heats []
    fun = @(b) myfun(t,b,g,M);                                              % Set given variables in myfun
    b   = fzero(fun,0.5);
    
    % Set solution variable based on input deg/rad
    if (strcmp(degrad,'deg'))                                               % If the output is in degrees
        sol = b*(180/pi);
    elseif (strcmp(degrad,'rad'))                                           % If the output is in radians
        sol  = b;
    end
    
end
%}

function [z0] = TM_Equations(theta,z,k)


% State-space form of the Taylor-Maccoll equation used for integration
% INPUTS
% theta : Integration angle [rad]
% z     : Angular and radial velocities [m/s]
% k     : Ratio of specific heats
% 
% OUTPUTS
% z0  : Solution array

% Initialize solution vector
z0 = zeros(2,1);                                                            % z0(1) = dVr/dTheta
                                                                            % z0(2) = d2Vr/dTheta2
% Define term used often in equation below
A = (k-1)/2;

% Numerator and denominator for z0 calculation below
num = (-2*A*z(1)) - (A*z(2)*cot(theta)) + (2*A*z(1)^3) + ...                % Numerator of second state-space equation
      (A*z(1)^2*z(2)*cot(theta)) + (2*A*z(1)*z(2)^2) + ...
      (A*z(2)^3*cot(theta)) + (z(1)*z(2)^2);
den = A*(1-z(1)^2-z(2)^2) - z(2)^2;                                         % Denominator of second state-space equation

% State-space representation
z0(1) = z(2);
z0(2) = num/den;

end

