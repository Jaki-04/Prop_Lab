function P = Intake_parameters()
    
    P.delta1 = deg2rad(11.110275451848205);
    P.delta2 = deg2rad(13.911672126164792)+P.delta1;
    P.delta3 = deg2rad(17.007645725250629)+P.delta2;
    P.alpha1 = deg2rad(24.6);
    P.alpha2 = deg2rad(31.6);
    P.alpha3 = deg2rad(44);
    P.Mf = 0.717801780013885;
    P.ptot_ratio = 0.763625294662585;
    P.p_ratio = 41.320623376265132;
end