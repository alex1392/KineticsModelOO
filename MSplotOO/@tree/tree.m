classdef tree
%% TREE  A class implementing a tree data structure.
%
% This class implements a simple tree data structure. Each node can only
% have one parent, and store any kind of data. The root of the tree is a
% privilieged node that has no parents and no siblings.
%
% Nodes are mainly accessed through their index. The index of a node is
% returned when it is added to the tree, and actually corresponds to the
% order of addition.
%
% Basic methods to tarverse and manipulate trees are implemented. Most of
% them take advantage of the ability to create _coordinated_ trees: If
% a tree is duplicated and only the new tree data content is modified
% (i.e., no nodes are added or deleted), then the iteration order and the
% node indices will be the same for the two trees.
%
% Internally, the class simply manage an array referencing the node parent
% indices, and a cell array containing the node data. 

% Jean-Yves Tinevez <tinevez@pasteur.fr> March 2012
    
    properties
        % Hold the data at each node
        Node = { [] };
        
        % Index of the parent node. The root of the tree as a parent index
        % equal to 0.
        Parent = [ 0 ]; %#ok<NBRAK>
        
    end
    
    methods
        
        % CONSTRUCTOR
        
        function [T, root_ID] = tree(content, val)
            %% TREE  Construct a new tree
            %
            % t = TREE(another_tree) is the copy-constructor for this
            % class. It returns a new tree where the node order and content
            % is duplicated from the tree argument.
            % 
            % t = TREE(another_tree, 'clear') generate a new copy of the
            % tree, but does not copy the node content. The empty array is
            % put at each node.
            %
            % t = TREE(another_tree, val) generate a new copy of the
            % tree, and set the value of each node of the new tree to be
            % 'val'.
            %
            % t = TREE(root_content) where 'root_content' is not a tree,
            % initialize a new tree with only the root node, and set its
            % content to be 'root_content'.
           
            if nargin < 1
                root_ID = 1;
                return
            end
            
            if isa(content, 'tree')
                % Copy constructor
                T.Parent = content.Parent;
                if nargin > 1 
                    if strcmpi(val, 'clear')
                        T.Node = cell(numel(T.Parent), 1);
                    else
                        cellval = cell(numel(T.Parent), 1);
                        for i = 1 : numel(T.Parent)
                            cellval{i} = val;
                        end
                        T.Node = cellval;
                    end
                else
                    T.Node = content.Node;
                end
                
            else
                % New object with only root content
                
                T.Node = { content };
                root_ID = 1;
            end
            
        end
        
        
        % METHODS
        
        function T = insertchild(T, parent, data, n)
          % n: insert the child as the n-th child
          childrenID = T.getchildren(parent);
          assert(n > 0 && n <= numel(childrenID) + 1,'Cannot insert n-th child for parentID = %d', parent);
          T = T.addnode(parent, data);
          
          childrenID = T.getchildren(parent);
          Nodes = T.Node(childrenID);
          Nodes = insert(Nodes(1:end-1), n, Nodes(end));
          T.Node(childrenID) = Nodes;
          
          for i = numel(childrenID) - 1:-1:n
            T.Parent(T.Parent == childrenID(i)) = childrenID(i+1);
          end
        end
        
        function [T, ID] = addnode(T, parent, data)
            %% ADDNODE attach a new node to a parent node
            % 
            % tree = tree.ADDNODE(parent_index, data) create a new node
            % with content 'data', and attach it as a child of the node
            % with index 'parent_index'. Return the modified tree.
            % 
            % [ tree ID ] = tree.ADDNODE(...) returns the modified tree and
            % the index of the newly created node.
            
            if parent < 0 || parent > numel(T.Parent)
                error('MATLAB:tree:addnode', ...
                    'Cannot add to unknown parent with index %d.\n', parent)
            end
            
            if parent == 0
                % Replace the whole tree by overiding the root.
                T.Node = { data };
                T.Parent = 0;
                ID = 1;
                return
            end
            
            % Expand the cell by
            T.Node{ end + 1, 1 } = data;
            
            T.Parent = [
                T.Parent
                parent ];
            
            ID = numel(T.Node);
        
        end
        
        function flag = isleaf(T, Idx)
           %% ISLEAF  Return true if given ID matches a leaf node.
           % A leaf node is a node that has no children.
           n = length(Idx);
           flag = true(1,n);
           for i = 1:n
               ID = Idx(i);
               assert( ID >= 1 && ID <= numel(T.Parent), 'MATLAB:tree:isleaf No node with ID %d.', ID)
               parent = T.Parent;
               flag(i) = ~any( parent == ID );
           end
           
        end
        
        function IDs = findleaves(T)
           %% FINDLEAVES  Return the IDs of all the leaves of the tree.
           parents = T.Parent;
           IDs = (1 : numel(parents)); % All IDs
           IDs = setdiff(IDs, parents); % Remove those which are marked as parent
           
        end
        
        function data = get(T, ID, toend)
            %% GET  Return the data of the given node ID.
            if strcmpi(ID,'all')
                ID = 1:T.nnodes;
            elseif strcmpi(ID,'end')
                ID = T.nnodes;
            end
            if nargin >= 3 && toend
                ID = ID:T.nnodes;
            end
            n = length(ID);
            data = [];
            for i = 1:n
                temp = T.Node{ID(i)};
                if ~isempty(temp)
                    data = [data, temp];
                end
            end
        end

        function T = set(T, ID, data)
            %% SET  Set the data of given node ID and return the modifed tree.
            T.Node{ID} = data;
        end
        
        function T = setsub(T,ID,sT)
          cID = [ID,T.getoffspring(ID)];
          sID = cID - ID + 1;
          for i = 1:numel(cID)
            T = T.set(cID(i),sT.get(sID(i)));
          end
        end
        
        function cID = getoffspring(T,ID)
          cID = [];
          cidx = T.getchildren(ID);
          if ~isempty(cidx)
            cID = [cID,cidx,T.getoffspring(cidx)];
          end
        end

        function cID = getchildren(T, ID)
            %% GETCHILDREN  Return the list of ID of the children of the given node ID.
            % The list is returned as a line vector.
            cID = [];
            n = length(ID);
            for i = 1:n
                Idx = ID(i);
                parent = T.Parent;
                temp = find(parent == Idx);
                cID = [cID,temp'];
            end
        end
        
        function ID = getparent(T,ID)
        %% GETPARENT  Return the ID of the parent of the given node.
            if ID < 1 || ID > numel(T.Parent)
                error('MATLAB:tree:getparent', ...
                    'No node with ID %d.', ID)
            end
            ID = T.Parent(ID);
        end
        
        function IDs = getancestors(T,ID)
          IDs = [];
          PID = T.getparent(ID);
          if PID
            IDs = [IDs,PID,T.getancestors(PID)];
          end
        end
        
        function IDs = getsiblings(T, ID)
            %% GETSIBLINGS  Return the list of ID of the sliblings of the 
            % given node ID, including itself.
            % The list is returned as a column vector.
            if ID < 1 || ID > numel(T.Parent)
                error('MATLAB:tree:getsiblings', ...
                    'No node with ID %d.', ID)
            end
            
            if ID == 1 % Special case: the root
                IDs = 1;
                return
            end
            
            parent = T.Parent(ID);
            IDs = T.getchildren(parent);
        end
        
        function n = nnodes(T)
            n = numel(T.Parent);
        end
        
        function n = n(T)
            n = numel(T.Parent);
        end
        
        function leaves = subleaves(T, ID)
            leaves = [];
            cidx = T.getchildren(ID);
            if ~isempty(cidx)
              leaves = [leaves, T.subleaves(cidx)];
            else
              leaves = ID;
            end
        end
        
    end
    
    % STATIC METHODS
    
    methods (Static)
        
        hl = decorateplots(ha)
        
        function [lineage, duration] = example
            
            lineage_AB = tree('AB');
            [lineage_AB, id_ABa] = lineage_AB.addnode(1, 'AB.a');
            [lineage_AB, id_ABp] = lineage_AB.addnode(1, 'AB.p');
            
            [lineage_AB, id_ABal] = lineage_AB.addnode(id_ABa, 'AB.al');
            [lineage_AB, id_ABar] = lineage_AB.addnode(id_ABa, 'AB.ar');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABal, 'AB.ala');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABal, 'AB.alp');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABar, 'AB.ara');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABar, 'AB.arp');
            
            [lineage_AB, id_ABpl] = lineage_AB.addnode(id_ABp, 'AB.pl');
            [lineage_AB, id_ABpr] = lineage_AB.addnode(id_ABp, 'AB.pr');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABpl, 'AB.pla');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABpl, 'AB.plp');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABpr, 'AB.pra');
            [lineage_AB, ~] = lineage_AB.addnode(id_ABpr, 'AB.prp');
            
            lineage_P1 = tree('P1');
            [lineage_P1, id_P2] = lineage_P1.addnode(1, 'P2');
            [lineage_P1, id_EMS] = lineage_P1.addnode(1, 'EMS');
            [lineage_P1, id_P3] = lineage_P1.addnode(id_P2, 'P3');
            
            [lineage_P1, id_C] = lineage_P1.addnode(id_P2, 'C');
            [lineage_P1, id_Ca] = lineage_P1.addnode(id_C, 'C.a');
            [lineage_P1, ~] = lineage_P1.addnode(id_Ca, 'C.aa');
            [lineage_P1, ~] = lineage_P1.addnode(id_Ca, 'C.ap');
            [lineage_P1, id_Cp] = lineage_P1.addnode(id_C, 'C.p');
            [lineage_P1, ~] = lineage_P1.addnode(id_Cp, 'C.pa');
            [lineage_P1, ~] = lineage_P1.addnode(id_Cp, 'C.pp');
            
            [lineage_P1, id_MS] = lineage_P1.addnode(id_EMS, 'MS');
            [lineage_P1, ~] = lineage_P1.addnode(id_MS, 'MS.a');
            [lineage_P1, ~] = lineage_P1.addnode(id_MS, 'MS.p');
            
            [lineage_P1, id_E] = lineage_P1.addnode(id_EMS, 'E');
            [lineage_P1, id_Ea] = lineage_P1.addnode(id_E, 'E.a');
            [lineage_P1, ~] = lineage_P1.addnode(id_Ea, 'E.al'); %#ok<*NASGU>
            [lineage_P1, ~] = lineage_P1.addnode(id_Ea, 'E.ar');
            [lineage_P1, id_Ep] = lineage_P1.addnode(id_E, 'E.p');
            [lineage_P1, ~] = lineage_P1.addnode(id_Ep, 'E.pl');
            [lineage_P1, ~] = lineage_P1.addnode(id_Ep, 'E.pr');
            
            [lineage_P1, id_P4] = lineage_P1.addnode(id_P3, 'P4');
            [lineage_P1, ~] = lineage_P1.addnode(id_P4, 'Z2');
            [lineage_P1, ~] = lineage_P1.addnode(id_P4, 'Z3');
            
            
            [lineage_P1, ~] = lineage_P1.addnode(id_P3, 'D');
            
            lineage = tree('Zygote');
            lineage = lineage.graft(1, lineage_AB);
            lineage = lineage.graft(1, lineage_P1);

            
            duration = tree(lineage, 'clear');
            iterator = duration.DFS;
            for i = iterator
               duration = duration.set(i, round(20*rand)); 
            end
            
        end
        
    end
    
end

