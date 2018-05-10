classdef Const
    
    properties (Constant)
        Tol = 1e-5;
        MTol = 1e-3;
        i_ini = 0.05;
        i_ang = [0, 1, sqrt(2)];
        A = [1000,1]; %[ incompatible, compatible ]
        vn = 6;        
        PlotL = 2;
        deffine = [2,5];
        MProp = MProp(1:6);
        varcolor = [70 130 180 ;...         % blue
            176 196 222;...                 % light blue
            139 125 107;...                 % brown
            205 175 149;...                 % light brown
            176 48 96;...                   % red
            255 192 203]/255;               % pink
    end
end