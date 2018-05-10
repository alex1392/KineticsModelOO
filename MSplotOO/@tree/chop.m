function T = chop(T,ID)
%% CHOP  Remove the target ID and all subIDs from the given tree.

    iterator = T.DFS(ID);
    
    % Build new parent array
    np = T.Parent;
    
    % Remove unwanted IDs
    np ( iterator ) = [];
    
    % Shift parent value: if a parent were after some IDs we removed, we
    % need to shift its value by an amount equal to the number of parent we
    % removed, and that were BEFORE the target parent
    for i = 1 : numel(np)
        np(i) = np(i) - sum(np(i) > iterator);
    end
    
    T.Parent = np;
    T.Node(iterator) = [];

end