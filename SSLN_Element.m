classdef SSLN_Element < handle
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Element class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        element_nodes % vector of two Node objects
        E % modulus of elasticity
        v % 
        A % area 
        Ayy % cross section shear area
        Azz % cross section shear area
        Izz % inertia about zz
        Iyy % inertia about yy
        webdir % 3x1 vector
        J % torsional section
        Cw % torsional constant
        node_i
        node_j
        element_length
        GAMMA
        k_local
        k_global
        
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %  Replace XYZ by your initials before proceeding
        function self = SSLN_Element(element_nodes,E,v,A,Ayy,Azz,Izz,Iyy,webdir,J)
            self.element_nodes = element_nodes;
            self.E = E;
            self.v = v;
            self.A = A;
            self.Ayy = Ayy;
            self.Azz = Azz;
            self.Izz = Izz;
            self.Iyy = Iyy;
            self.webdir = webdir;
            self.J = J;
            ComputeLength(self);
            ComputeTransformationMatrix(self)
            ComputeElasticStiffnessMatrix(self)
        end
        
        function  length = GetNodeLength(self)
            length = self.element_length;
        end
        function  [GAMMA] = GetNodegamma(self)
            
            GAMMA = self.GAMMA;
        end
            
            
    end
    
    % Private methods go here
    methods (Access = private)
        %% Compute the element's length
        function ComputeLength(self)
            self.node_i = self.element_nodes(1); %N1
            self.node_j = self.element_nodes(2);
            i_coords = GetNodeCoord(self.node_i);
            j_coords = GetNodeCoord(self.node_j);
            self.element_length = sqrt((j_coords(1) - i_coords(1))^2 + (j_coords(2) - i_coords(2))^2 + (j_coords(3) - i_coords(3))^2);
            disp(self.element_length);
        end
        
        %% Compute the element's geometric transformation matrix
        function ComputeTransformationMatrix(self)
            i_coords = GetNodeCoord(self.node_i);
            j_coords = GetNodeCoord(self.node_j);
            lamda_x_prime = (j_coords(1) - i_coords(1))/self.element_length;
            u_x_prime = (j_coords(2) - i_coords(2))/self.element_length;
            v_x_prime = (j_coords(3) - i_coords(3))/self.element_length;
            x_prime = [lamda_x_prime; u_x_prime; v_x_prime];
            z_prime = cross(x_prime, (self.webdir/norm(self.webdir)).');
            gamma = [x_prime.'; self.webdir/norm(self.webdir); z_prime.'];
            zero = zeros(3);
            self.GAMMA = [gamma, zero, zero, zero; ...
                          zero, gamma, zero, zero; ...
                          zero, zero, gamma, zero; ...
                          zero, zero, zero, gamma];
            disp(self.GAMMA)
        end
        
        %% Compute the element's elastic stiffness matrix in local and global coordinates
        function ComputeElasticStiffnessMatrix(self)
            a = self.A / self.element_length;
            b = 12 * self.Izz / (self.element_length^3);
            c = 6 * self.Izz / (self.element_length^2);
            d = 12 * self.Iyy / (self.element_length^3);
            e = 6 * self.Iyy / (self.element_length^2);
            f = self.J / (2 * (1 + self.v) * self.element_length);
            g = 4 * self.Iyy / self.element_length;
            h = 2 * self.Iyy / self.element_length;
            i = 4 * self.Izz / self.element_length;
            j = 2 * self.Izz / self.element_length;
            
            self.k_local = self.E * [a, 0, 0, 0, 0, 0, -a, 0, 0, 0, 0, 0; ...
                                     0, b, 0, 0, 0, c, 0, -b, 0, 0, 0, c; ...
                                     0, 0, d, 0, -e, 0, 0, 0, -d, 0, -e, 0; ...
                                     0, 0, 0, f, 0, 0, 0, 0, 0, -f, 0, 0; ...
                                     0, 0, -e, 0, g, 0, 0, 0, e, 0, h, 0; ...
                                     0, c, 0, 0, 0, i, 0, c, 0, 0, 0, j; ...
                                     -a, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0; ...
                                     0, -b, 0, 0, 0, c, 0, b, 0, 0, 0, -c; ...
                                     0, 0, -d, 0, e, 0, 0, 0, d, 0, e, 0; ...
                                     0, 0, 0, -f, 0, 0, 0, 0, 0, f, 0, 0; ...
                                     0, 0, -e, 0, h, 0, 0, 0, e, 0, g, 0; ...
                                     0, c, 0, 0, 0, j, 0, -c, 0, 0, 0, i];
                                 
            self.k_global = self.GAMMA.' * self.k_local * self.GAMMA;
            disp(self.k_local)
            disp(self.k_global)
        end
    end
end
