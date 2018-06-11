function TA = graft(TA, ID, TB)
%% GRAFT   Graft another tree at the specified node of this tree.

    
    nNodes = numel(TA.Parent);

    otParents = TB.Parent;
    % Shift other parent indices
    otParents = otParents + nNodes;
    % Make the other root a child of the target node
    otParents(1) = ID;
    
    % Concatenate
    newParents = [ TA.Parent ; otParents ];
    newNodes   = vertcat( TA.Node, TB.Node );
    
    % Edit
    TA.Node = newNodes;
    TA.Parent = newParents;
    

end