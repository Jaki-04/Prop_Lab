function Parameters = Geometry()
    Parameters.delta1 = deg2rad(11.8);
    Parameters.delta2 = deg2rad(14.4+11.8);
    Parameters.delta3 = deg2rad(17.7+14.4+11.8);
    Parameters.l1 = 0.428231085732875;          % Prima rampa
    Parameters.l2 = 0.264764603971497;          % Seconda rampa
    Parameters.l3_shock = 0.159533184247144;    % Terza rampa (fino all'urto normale)
    Parameters.l3 = 0.166738695363947;          % Terza rampa (fino alla sezione di ingresso supersonica)
    Parameters.L_sup = 0.700536761338328;       % Distanza orizzontale dalla punta del cono al cowl lip
    Parameters.L_sub = 0.914613896582030;       % Distanza orizzontale tra il cowl lip e il compressore
    Parameters.ds = 0.01;                       % Distanza tra l'onda normale e la sezione di ingresso (parallela a rampa 3)

    Parameters.Ain_sup = 0.667777653078562;         % Area di ingresso in supersonico
    Parameters.Hin_sup = 0.229589275998660;         % Altezza all'ingresso in supersonico
    Parameters.rin = 0.380198703067846;             % Raggio minimo all'ingresso in supersonico
    Parameters.Rin = 0.545629511118563;             % Raggio massimo all'ingresso in supersonico
    Parameters.A1_sup =  1.514674633685526;         % Area dopo la diffusione subsonica

    Parameters.Ain_sub = 0.452820254995432;         % Area di ingresso in subsonico
    Parameters.Hin_sub = 0.152624354077758;         % Altezza all'ingresso in subsonico
    Parameters.A1_sub = 0.602370359775345;          % Area dopo la diffusione

    Parameters.Rc_h1 = 0.253280069296420;       % Raggio hub del compressore
    Parameters.Rc_t = 0.506560138592840;        % Raggio tip del compressore
    Parameters.Rc_he = 0.474852161629951;       %  Raggio finale dell'hub
end