clear 
close all
clc

%% Dimensionamento post combustore

supCruise = Supercruise_parameters();

m_a = supCruise.m_a;
f = supCruise.f;

% da presa.m da sistemare con mach 0,25
Ttot1 = 747.4425;
v2 = 267.4095;
A2 = 0.8317; % da sistemare perchè dovrebbe esserci pi_AB = 0.95 se acceso
r_t5 = 0.45; %numero a caso da calcolare poi
L_H_ratio = 3.2; %arbitrario
B = 0.28; %arbitrario

% Parametri termodinamici
cp_a = 1005;                % Calore specifico aria [J/kgK] in air parameters...
                            % è 1000, secondo me deve essere 1005
cp_GC = 1150;               % Calore specifico gas combusti [J/kgK]
H_f = 43000000;             % Potere calorifico [J/kg]
eta_AB = 0.92;              % Rendimento post-combustore

% Temperatura all'uscita
T_t7 = (cp_a * Ttot1 + f * H_f * eta_AB) / ((1 + f) * cp_GC);

% dimensioni post bruciatore
r_ram = sqrt((A2 / pi) + r_t5^2);
A_ram = pi * r_ram^2;
H_ram = 2 * r_ram;
L_ram = L_H_ratio * H_ram;

% dimensioni stabilizzatori di fiamma (ne metto 2)
d_fh = B * H_ram / 2;