classdef Const < handle
  
  properties (Constant)
    FerroConst = Const
  end
  
  properties
    Tol
    MTol
    i_ini
    i_ang
    A
    vn
    PlotL
    deffine
    MProp
    varcolor
  end
  
  methods (Access = private)
    function obj = Const %must be the constructor!!
      obj.Tol = 1e-5;
      obj.MTol = 1e-3;
      obj.i_ini = 0.05;
      obj.i_ang = [0, 1, sqrt(2)];
      obj.A = [1000,1]; %[ incompatible, compatible ]
      obj.vn = 6;
      obj.PlotL = 2;
      obj.deffine = [2,5];
      obj.MProp = MProp(1:6);
      obj.varcolor = [70 130 180 ;...         % blue
        176 196 222;...                 % light blue
        139 125 107;...                 % brown
        205 175 149;...                 % light brown
        176 48 96;...                   % red
        255 192 203]/255;               % pink
    end
  end
end