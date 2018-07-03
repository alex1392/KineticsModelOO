function multirun
% var = {[5 3 5 4],[5 3 6 1],[5 6 5 3],[5 6 6 5],[5 5 5 6],[5 3 4 6]};
sigma3 = -1.78;%,[0;0;-0.36;0;0;0],[0;0;-0.72;0;0;0],[0;0;-1.07;0;0;0],[0;0;0;0;0;0]};
EVec3 = 1;
startpath = pwd;
immobility = 0.01;
for i = 1:5
  A = [0,1]*10^i;
  cd(startpath)
  try
    KineticsModel(sigma3,EVec3,immobility,A)
  catch
    continue
  end
end