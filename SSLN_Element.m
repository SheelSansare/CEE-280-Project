classdef SSLN_Element < handle

% Element class for a 3-dimensional framed structure
    
    % Private properties 
    properties (Access = private)
        ele_dofs
        element_length
        GAMMA
        k_local
        k_global
        local_FEF
        global_FEF
    end
    
    % Public methods
    methods (Access = public)
        %% Constructor
        function self = SSLN_Element(element_nodes,E,v,A,Ayy,Azz,Izz,Iyy,webdir,J, w)
            ComputeLength(self, element_nodes);
            ComputeTransformationMatrix(self, webdir, element_nodes);
            ComputeElasticStiffnessMatrix(self, E,v,A,Ayy,Azz,Izz,Iyy,J);
            ComputeFixedEndForces(self, w);
            RetrieveDOFS(self, element_nodes);
        end
        
        %% Getter for Global Stiffness Matrix
        function k_global = GetGlobalStiffness(self)
            k_global = self.k_global;
        end
        
        %% Getter for Element DOFS
        function ele_dofs = GetDOFS(self)
            ele_dofs = self.ele_dofs;
        end
        
        %% Getter For Fixed End Forces
        function global_FEF = GetGlobalFEF(self)
            global_FEF = self.global_FEF;
        end
        
        %% Compute and Return Forces (Force Recovery)
        function f_local = ComputeForces(self, u)
            f_global = self.k_global * u;
            f_local = (self.GAMMA * f_global) + self.local_FEF;
        end
    end

    % Private methods 
    methods (Access = private)
        %% Assemble Element Degrees of Freedom 
        function RetrieveDOFS(self, element_nodes)
            node_i = element_nodes(1);
            node_j = element_nodes(2);
            i_dofs = GetNodeDofs(node_i);
            j_dofs = GetNodeDofs(node_j);
            self.ele_dofs = [i_dofs; j_dofs];
        end
        %% Compute Element's Length
        function ComputeLength(self, element_nodes)
            node_i = element_nodes(1);
            node_j = element_nodes(2);
            i_coords = GetNodeCoord(node_i);
            j_coords = GetNodeCoord(node_j);
            self.element_length = sqrt((j_coords(1) - i_coords(1))^2 + ...
                (j_coords(2) - i_coords(2))^2 + (j_coords(3) - i_coords(3))^2);
        end
        
        %% Compute Element's Geometric Transformation Matrix
        function ComputeTransformationMatrix(self, webdir, element_nodes)
            node_i = element_nodes(1);
            node_j = element_nodes(2);
            i_coords = GetNodeCoord(node_i);
            j_coords = GetNodeCoord(node_j);
            lamda_x_prime = (j_coords(1) - i_coords(1))/self.element_length;
            u_x_prime = (j_coords(2) - i_coords(2))/self.element_length;
            v_x_prime = (j_coords(3) - i_coords(3))/self.element_length;
            x_prime = [lamda_x_prime; u_x_prime; v_x_prime];
            z_prime = cross(x_prime, (webdir/norm(webdir)).');
            gamma = [x_prime.'; webdir/norm(webdir); z_prime.'];
            zero = zeros(3);
            self.GAMMA = [gamma, zero, zero, zero; ...
                          zero, gamma, zero, zero; ...
                          zero, zero, gamma, zero; ...
                          zero, zero, zero, gamma];
        end
        
        %% Compute Element's Elastic Stiffness Matrix in Local and Global Coordinates
        function ComputeElasticStiffnessMatrix(self, E,v,A,Ayy,Azz,Izz,Iyy,J)
            G = E/2/(1+v);
            phi_y = 12*E*Izz/G/Ayy/self.element_length^2;
            phi_z = 12*E*Iyy/G/Azz/self.element_length^2;
            a = A / self.element_length;
            b = 12 * Izz / (self.element_length^3) / (1+phi_y);
            c = 6 * Izz / (self.element_length^2) / (1+phi_y);
            d = 12 * Iyy / (self.element_length^3) / (1+phi_z);
            e = 6 * Iyy / (self.element_length^2) / (1+phi_z);
            f = G * J / self.element_length / E;
            g = (4 + phi_z) * Iyy / self.element_length /(1 + phi_z);
            h = (2 - phi_z) * Iyy / self.element_length / (1 + phi_z);
            i = (4 + phi_y) * Izz / self.element_length /(1 + phi_y);
            j = (2 - phi_y) * Izz / self.element_length / (1 + phi_y);

            self.k_local = E * [a, 0, 0, 0, 0, 0, -a, 0, 0, 0, 0, 0; ...
                                     0, b, 0, 0, 0, c, 0, -b, 0, 0, 0, c; ...
                                     0, 0, d, 0, -e, 0, 0, 0, -d, 0, -e, 0; ...
                                     0, 0, 0, f, 0, 0, 0, 0, 0, -f, 0, 0; ...
                                     0, 0, -e, 0, g, 0, 0, 0, e, 0, h, 0; ...
                                     0, c, 0, 0, 0, i, 0, -c, 0, 0, 0, j; ...
                                    -a, 0, 0, 0, 0, 0, a, 0, 0, 0, 0, 0; ...
                                     0, -b, 0, 0, 0, -c, 0, b, 0, 0, 0, -c; ...
                                     0, 0, -d, 0, e, 0, 0, 0, d, 0, e, 0; ...
                                     0, 0, 0, -f, 0, 0, 0, 0, 0, f, 0, 0; ...
                                     0, 0, -e, 0, h, 0, 0, 0, e, 0, g, 0; ...
                                     0, c, 0, 0, 0, j, 0, -c, 0, 0, 0, i];
                                 
            self.k_global = self.GAMMA.' * self.k_local * self.GAMMA;
        end
        
        %% Compute Fixed End Forces
        function ComputeFixedEndForces(self, w)
            self.local_FEF = [-w(1)*self.element_length/2;
                              -w(2)*self.element_length/2;
                              -w(3)*self.element_length/2;
                               0;
                               w(3)*self.element_length^2/12;
                              -w(2)*self.element_length^2/12;
                              -w(1)*self.element_length/2;
                              -w(2)*self.element_length/2;
                              -w(3)*self.element_length/2;
                               0;
                              -w(3)*self.element_length^2/12;
                               w(2)*self.element_length^2/12];
                     
            self.global_FEF = self.GAMMA.' * self.local_FEF;
        end  
    end
end
