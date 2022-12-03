classdef SSLN_Node < handle

% Node class for a 3-dimensional framed structure
    
    % Private properties 
    properties (Access = private)
        node_coord
        
        dofs
    end
    
    % Public methods 
    methods (Access = public)
        %% Constructor
        function self = SSLN_Node(node_coord, node_num)
            self.node_coord = node_coord;
            AssignDOF(self, node_num);
        end
        
        %% Getter Node Coordinates
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
        
        %% Getter for Node Degrees of Freedom
        function node_dofs = GetNodeDofs(self)
            node_dofs = self.dofs;
        end
    end
    
    % Private methods
    methods (Access = private)
        %% Assign Degrees of Freedom to Node
        function AssignDOF(self, node_num)
            %creates empty vector
            self.dofs = NaN(6,1);
            index = 5;
            for i = 1:6
                self.dofs(i, 1) = node_num*6 - index;
                index = index - 1;
            end
        end
    end
end
