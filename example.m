addpath(genpath(pwd)) % Add ./MSPlotOO into path

var = [5 3 5 4]; %define variants
dof = [0.5 0.5]; %define degree of freedoms
type = 'COA'; %define cut type
fine = 2; %define fineness
obj = MStruct(var,dof,type,'fine',fine); %Construct microstructure object
[obj, MST, IFT] = obj.getTree; %Construct tree diagram for the object

hMS = MST.draw; %draw the microstructure
