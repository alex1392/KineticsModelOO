function multirun
% var = {[5 3 5 4],[5 3 6 1],[5 6 5 3],[5 6 6 5],[5 5 5 6],[5 3 4 6]};
sigma = {[0;0;-1.78;0;0;0],[0;0;-0.36;0;0;0]};%[0;0;-0.72;0;0;0],[0;0;-1.07;0;0;0],
for i = 1:numel(sigma)
    try 
      KineticsModel(sigma{i})
      cd ..
    catch
      cd ..
      continue
    end
end