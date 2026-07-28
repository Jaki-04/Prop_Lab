function Parameters = Geometry()

% Spina supersonica
    Spina = Intake_parameters();
    Parameters.delta1 = Spina.delta1;           % Angolo della prima rampa
    Parameters.delta2 = Spina.delta2;           % Angolo della seconda rampa
    Parameters.delta3 = Spina.delta3;           % Angolo della terza rampa
    Parameters.L_sup = 0.754408660564071;       % Distanza orizzontale dalla punta del cono al cowl lip
    Parameters.ds = 0.01;                       % Distanza tra l'onda normale e la sezione di ingresso (parallela a rampa 3)

    Parameters.Ain_sup = 0.646624148653795;         % Area di ingresso in supersonico
    Parameters.Hin_sup = 0.223388871860380;         % Altezza all'ingresso in supersonico
    Parameters.rin = 0.377725310475119;             % Raggio minimo all'ingresso in supersonico
    Parameters.Rin = 0.543658367915731;             % Raggio massimo all'ingresso in supersonico
    Parameters.A1_sup =  1.450055123493309;         % Area dopo la diffusione subsonica del flusso supersonico
 
% Diffusore subsonico
    Parameters.Ain_sub = 0.471664627794136;         % Area di ingresso in subsonico
    Parameters.Hin_sub = 0.162306945976957;         % Altezza all'ingresso in subsonico
    Parameters.r_sub = 0.381351421938773;           % Raggio di sommità della spina
    Parameters.A1_sub = 0.627438345355209;          % Area dopo la diffusione
    Parameters.phi_sup = deg2rad(-2.471455796);     % Angolo superiore di divergenza (positivo se antiorario)
    Parameters.phi_inf = deg2rad(-10);              % Angolo inferiore di divergenza (positivo se antiorario)
    Parameters.WLratio = 4.360593112483700;         % Rapporto L/W del diffusore
    Parameters.L_sub = 0.760218794189371;           % Distanza orizzontale tra il cowl lip e il compressore

% Compressore
    Parameters.Rc_h1 = 0.2422;                  % Raggio iniziale dell'hub del compressore
    Parameters.Rc_t = 0.4843;                   % Raggio tip del compressore
    Parameters.Rc_he = 0.4537;       % Raggio finale dell'hub del compressore
    Parameters.L_comp = 1.033595957996741;      % Lunghezza del compressore

% Camera di combustione
    Parameters.L_diff = 0.090322521873416;          % Lunghezza del diffusore
    Parameters.phi_diff = deg2rad(3);               % Angolo di divergenza
    Parameters.R_diff = 0.522336533245617;          % Raggio massimo del diffusore
    Parameters.r_diff = 0.474592088960632;          % Raggio minimo del diffusore
    Parameters.WLratio_diff = 2.570121914518049;    % Rapporto L/W del diffusore

    Parameters.rm = 0.498464311103125;          % Raggio medio della camera
    Parameters.href = 0.145861686586823;        % Spessore della camera
    Parameters.Aref = 0.456830580994646;        % Area di una sezione della camera
    Parameters.Lcc = 0.670414045111986;         % Lunghezza della camera
    
end