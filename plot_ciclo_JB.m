%% ==========================================================================
%  DIAGRAMMA h-s DI UN CICLO JOULE-BRAYTON (turbogetto semplice)
%  Fasi:
%   0s      : ambiente (stato statico)
%   0s -> 0 : accelerazione/decelerazione isoentropica (retta verticale tratteggiata)
%   0 -> 2  : diffusore adiabatico NON isoentropico
%   2 -> 3  : compressore adiabatico NON isoentropico
%   3 -> 4  : combustione a pressione costante (isobara)
%   4 -> 5  : turbina adiabatica NON isoentropica
%   5 -> 7  : postcombustore spento (sola perdita di piccola entita')
%   7 -> 9  : ugello di scarico isoentropico
% ==========================================================================

clear; clc; close all;

%% --- Parametri termodinamici schematici -----------------------------------
cp    = 1.3;   % [kJ/(kg K)] calore specifico gas combusti (schematico, costante)
Href  = 80;    % [kJ/kg] offset per costruire correttamente l'isobara 3-4

% Entita' delle irreversibilita' (NON esagerate)
ds_diff = 0.06;  % aumento di entropia nel diffusore
ds_comp = 0.06;  % aumento di entropia nel compressore
ds_turb = 0.035;  % aumento di entropia in turbina
ds_duct = 0.020;  % aumento di entropia nel postcombustore spento

%% --- Punti principali del ciclo (s in kJ/(kg K), h in kJ/kg) --------------
P0s = [0.000,   0];                      % ambiente statico
P0  = [0.000,  1*(216.65*(1+0.2*3.5^2))];                      % ristagno ambiente (stessa s di 0s)

P2  = [P0(1) + ds_diff,  248];            % uscita diffusore
P3  = [P2(1) + ds_comp, 500];            % uscita compressore

h4  = 1*(800);                                % uscita camera di combustione
s4  = cp*log((h4-P3(2))/P3(2)+1)+P3(1);     % s coerente con isobara p3=p4
P4  = [s4, h4];

P5  = [P4(1) + ds_turb, 500];            % uscita turbina (reale)
P7  = [P5(1) + ds_duct, 500];            % uscita postcombustore spento
P9  = [P7(1), 100];                      % uscita ugello (isoentropico -> s7=s9)

%% --- Curva isobara di combustione (3 -> 4), forma esponenziale reale ------
%s_comb = linspace(P3(1), P4(1), 60);
%h_comb = (P3(2)) .* exp( (s_comb - P3(1)) / cp );

%% --- Curve di pressione totale costante
function s = sval(state, amp1, amp2)
if nargin==2
    amp2=amp1;
end
    s=linspace(state(1)-amp1, state(1)+amp2,60);
end
function h = hs(sval, state)
    cp = 1.3;
    h=state(2)+state(2)*(exp(( sval-state(1) )/cp) - 1 );
end
s_0s = sval(P0s, 0.02,0.1);
s_0 = sval(P0, 0.05);
s_2 = sval(P2,0.05);
s_3 = sval(P3,0.05, 0.95);
s_4 = sval(P4,0.05);
s_5 = sval(P5,0.05);
s_7 = sval(P7,0.05);

s34 = sval(P3, 0, P4(1)-P3(1));
h34 = hs(s34, P3)-(s34-s34(1))*30;
P4=[s34(end), h34(end)];

h_0s = hs(s_0s, P0s);
h_0 = hs(s_0, P0);
h_2 = hs(s_2, P2);
h_3 = hs(s_3, P3);
h_4 = hs(s_4, P4);
h_5 = hs(s_5, P5);
h_7 = hs(s_7, P7);


%% --- Figura ----------------------------------------------------------------
figure('Color','w'); hold on; box on; grid on;

col_red  = [0.80 0.10 0.10];
col_gray = [0.55 0.55 0.55];

plot(s_0s,h_0s,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_0,h_0,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_2,h_2,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_3,h_3,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_4,h_4,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_5,h_5,'-', 'Color', col_gray, 'LineWidth', 1.2)
plot(s_7,h_7,'-', 'Color', col_gray, 'LineWidth', 1.2)

% ---- Trasformazioni reali (ROSSE) -----------------------------------------
plot([P0s(1) P0(1)], [P0s(2) P0(2)], '--', 'Color', col_red, 'LineWidth', 1.4); % 0s-0 tratteggiata
plot([P0s(1)  P2(1)], [P0s(2)  P2(2)], '-',  'Color', col_red, 'LineWidth', 1.6); % 0-2 diffusore
plot([P2(1)  P3(1)], [P2(2)  P3(2)], '-',  'Color', col_red, 'LineWidth', 1.6); % 2-3 compressore
plot(s34, h34, '-',  'Color', col_red, 'LineWidth', 1.6); % 3-4 combustione
plot([P4(1)  P5(1)], [P4(2)  P5(2)], '-',  'Color', col_red, 'LineWidth', 1.6); % 4-5 turbina
plot([P5(1)  P7(1)], [P5(2)  P7(2)], '-',  'Color', col_red, 'LineWidth', 1.6); % 5-7 postcombustore spento
plot([P7(1)  P9(1)], [P7(2)  P9(2)], '-',  'Color', col_red, 'LineWidth', 1.6); % 7-9 ugello isoentropico

% ---- Trasformazioni isoentropiche di confronto (GRIGIE) -------------------
% (una per ogni trasformazione realmente non isoentropica, tra le stesse h)
%plot([P0(1) P0(1)], [P0(2) P2(2)], '-', 'Color', col_gray, 'LineWidth', 1.2); % diffusore ideale
%plot([P2(1) P2(1)], [P2(2) P3(2)], '-', 'Color', col_gray, 'LineWidth', 1.2); % compressore ideale
%plot([P4(1) P4(1)], [P4(2) P5(2)], '-', 'Color', col_gray, 'LineWidth', 1.2); % turbina ideale
%plot([P5(1) P5(1)], [P5(2) P7(2)], '-', 'Color', col_gray, 'LineWidth', 1.2); % postcombustore ideale

% ---- Punti e etichette ------------------------------------------------------
pts   = [P0s; P0; P2; P3; P4; P5; P7; P9; P0s; P0; P2; P3; P4; P5; P7; P9];
names = {'0s','0','2','3','4','5','7','9', 'P_{0s}', 'P_0', 'P_2', 'P_3', 'P_4', 'P_5', 'P_7', 'P_9'};

plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 4);
for k = 1:(size(pts,1)/2)
    text(pts(k,1) + 0.003, pts(k,2), names{k}, 'FontSize', 11, 'FontWeight','bold');
end
for k = (size(pts,1)/2):size(pts,1)
    text(pts(k,1) + 0.008, pts(k,2), names{k}, 'FontSize', 8, 'Color', col_gray);
end

xlabel('s');
ylabel('h');
%title('Diagramma h-s del ciclo Joule-Brayton');
xlim([-0.1, s4 + 0.2]);
ylim([-50, P4(2)+60])
set(gca, 'XTick', [], 'YTick', [])
%legend({'trasformazioni reali (e riferimento 0s-0)'}, 'Location','southeast');