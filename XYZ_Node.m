classdef XYZ_Node < handle
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Node class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        % 3x1 vector containing the x, y, and z coordinates of the node
        node_coord
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %  Replace XYZ by your initials before proceeding
        %    Arguments
        %      node_coord:  3x1 vector containing the x, y, and z coordinates of the node
        function self = XYZ_Node(node_coord)
            self.node_coord = node_coord;
        end
        
        %% Get Node Coordinates
        %  Return "node_coord"
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        
    end
end
