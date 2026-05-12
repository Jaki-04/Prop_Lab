function Air = Air_parameters(varargin)


Air.cp_GC=1243;
Air.g_GC=1.33;
Air.R_GC=Air.cp_GC/3.5;

Air.cp=1004;
Air.g=1.4;
Air.R=287.01;



if ismember('25000', varargin)
    Air.p= 2549; 
    Air.T= 216.65;
    Air.rho = Air.p/(Air.R*Air.T);
elseif ismember ('12000', varargin)
    Air.p= 19267; 
    Air.T= 216.65;
    Air.rho = Air.p/(Air.R*Air.T);
end
