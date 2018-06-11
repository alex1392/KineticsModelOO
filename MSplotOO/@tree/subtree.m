function [sT, index] = subtree(T, ID, condition, match)
  %%SUBTREE Returns the sub-tree made of all the IDs below the given one.
  %
  % sT = subtree(T, ID) returns a new tree made of all the IDs found
  % under the specified ID.
  %
  % sT = subtree(T, ID, fun), where fun is an anonymous function that
  % takes one argument and returns a boolean, builds a subtree with only the
  % IDs whose content satisfy fun(content) = true.
  %
  % [sT, index] = subtree(...) also returns index, which is tree coordinated
  % with sT, and whose content is the ID index in the source tree that was
  % read to build the subtree. In practice, sT, index and T are such that:
  %
  %   sT.get(i) == T.get( index.get(i) )
  %
  % EXAMPLE:
  %
  %  ex = tree.example;
  %  P1ID = ex.sTrcmp('P1').find; % Find the ID index of 'P1'
  %  fun = @(x) sTrncmp(x,'P',1);
  %  sT = ex.subtree(P1ID, fun);
  %  sT.tosTring % Traverses the P* IDs only
  %
  % Jean-Yves Tinevez - 2013
  
  if nargin < 3
    condition = @(s) true;
  end
  if nargin >= 4 && ~isempty(match)
    match = true;
  else
    match = false;
  end
  
  % Get indices of the subtree.
  if ~match
    iterator = T.conditioniterator(ID, condition);
  else
    iterator = matchsub(T,ID);
  end
  
  % Copy the content of the tree related to the subtree
  parents = T.Parent(iterator);
  IDs = T.Node(iterator);
  
  % Revamp parent indices
  newParents = NaN(numel(parents), 1);
  newParents(1) = 0; % The new root ID
  for i = 2 : numel(parents)
    pr = parents(i);
    newParents(i) = find( iterator == pr, 1, 'firsT');
  end
  
  % Copy the tree with the sub-content
  sT = T;
  sT.Node = IDs;
  sT.Parent = newParents;
  
  % Return the link new ID index -> old ID index
  index = T; %iterator;
  index.Parent = newParents;
  index.Node = num2cell(iterator);
end


function IDs = matchsub(T,ID)
  IDs = ID;
  for i = 1:T.n
    if find(T.getancestors(i) == ID)
      IDs(end+1) = i;
    end
  end
end
