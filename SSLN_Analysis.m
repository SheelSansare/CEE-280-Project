classdef SSLN_Analysis < handle

% Analysis class for a 3-dimensional framed structure
    
    % Private properties 
    properties (Access = private)
        node_vector
        element_vector
        global_concen
        global_fef
        freeDOFS
        fixedDOFS
        dispDOFS
        fix_t
        AFLAG
        DEFL
        REACT
        ELE_FOR
    end
    
    % Public methods 
    methods (Access = public)
        %% Constructor
        function self = SSLN_Analysis(nnodes,coord,nele,A,Izz,Iyy,J,Ayy,Azz,E,v,webdir,w,ends)
            CreateNodes(self, nnodes, coord);
            CreateElements(self, nele, ends,E,v,A,Ayy,Azz,Izz,Iyy,webdir,J, w);
        end
        
        %% Run the analysis
        function RunAnalysis(self, fixity, nnodes, concen)
            K = CreateStiffnessMatrix(self, nnodes);
            ClassifyDOFS(self, fixity);
            [ff, fn, fs, feff, fefn, fefs] = CreateLoadVector(self, nnodes, concen);
            [kff, kfn, kfs, knf, knn, kns, ksf, ksn, kss] = ComputeStiffnessSubMatrices(self, K);
            CheckKffMatrix(self, kff);
            if self.AFLAG == 0
                self.DELF = 0;
                self.REACT = 0;
                self.ELE_FOR = 0;
            else
                [self.DEFL, self.REACT] = ComputeDisplacementsReactions(self, kff, ff, feff, kfn, ksf, ksn, nnodes);
                self.ELE_FOR = RecoverElementForces(self, self.DEFL);
                ComputeError(self, kff, kfn, ff, feff);
            end
        end
        
        %% Get MASTAN results
        function [AFLAG, DEFL, REACT, ELE_FOR] = GetMastan2Returns(self)
            AFLAG = self.AFLAG;
            DEFL = self.DEFL;
            REACT = self.REACT;
            ELE_FOR = self.ELE_FOR;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Create Nodes vector
        function CreateNodes(self, nnodes, coord)
            % creates empty vector
            self.node_vector = [];
            for i = 1:nnodes
                self.node_vector = [self.node_vector, SSLN_Node(coord(i,:).', i)];
            end
        end
        
        %% Create Elements vector
        function CreateElements(self, nele, ends,E,v,A,Ayy,Azz,Izz,Iyy,webdir,J, w)
            % creates empty vector
            self.element_vector = [];
            for i = 1:nele
                element_nodes = [self.node_vector(ends(i,1)),self.node_vector(ends(i,2))];
                self.element_vector = [self.element_vector, ...
                    SSLN_Element(element_nodes,E(i),v(i),A(i),Ayy(i),Azz(i),Izz(i),...
                    Iyy(i),webdir(i,:),J(i), w(i,:))];
            end    
        end
        
        %% Create Global Stiffness Matrix
        function K = CreateStiffnessMatrix(self, nnodes)
            %creates empty matrix
            K = zeros(6*nnodes, 6*nnodes);
            for i = 1:length(self.element_vector)
                dofs = GetDOFS(self.element_vector(i));
                K(dofs, dofs) = K(dofs, dofs) + GetGlobalStiffness(self.element_vector(i));
            end
            % removes zeros for matrix
            K = sparse(K);
        end
        
        %% Create Load Vectors
        function [ff, fn, fs, feff, fefn, fefs] = CreateLoadVector(self, nnodes, concen)
            % global structure concentrated loads
            self.global_concen = transpose(concen(1,:));
            for i = 2:nnodes
                self.global_concen(1:6*i, 1) = [self.global_concen; transpose(concen(i,:))];
            end
            ff = self.global_concen(self.freeDOFS);
            fn = self.global_concen(self.dispDOFS);
            fs = self.global_concen(self.fixedDOFS);
            
            % global structure FEF
                % creates empty vector
            self.global_fef = zeros(6*nnodes, 1);
            for i = 1:length(self.element_vector)
                element_fefs = GetGlobalFEF(self.element_vector(i));
                dofs = GetDOFS(self.element_vector(i));
                self.global_fef(dofs) = self.global_fef(dofs) + element_fefs;
            end
            feff = self.global_fef(self.freeDOFS);
            fefn = self.global_fef(self.dispDOFS);
            fefs = self.global_fef(self.fixedDOFS);
        end
        
        %% Seperate and Store Degrees of Freedom
        function ClassifyDOFS(self, fixity)
            self.fix_t = transpose(fixity);
            self.freeDOFS = find(isnan(self.fix_t));
            self.fixedDOFS = find(self.fix_t == 0);
            self.dispDOFS = find(self.fix_t ~= 0 & ~isnan(self.fix_t));
        end
        
        %% Compute Stiffness Submatrices
        function [kff, kfn, kfs, knf, knn, kns, ksf, ksn, kss] = ComputeStiffnessSubMatrices(self, K)
            kff = K(self.freeDOFS, self.freeDOFS);
            kfn = K(self.freeDOFS, self.dispDOFS);
            kfs = K(self.freeDOFS, self.fixedDOFS);
            knf = transpose(kfn);
            knn = K(self.dispDOFS, self.dispDOFS);
            kns = K(self.dispDOFS, self.fixedDOFS);
            ksf = transpose(kfs);
            ksn = transpose(kns);
            kss = K(self.fixedDOFS, self.fixedDOFS);
        end
        
        %% Check Kff matrix
        function CheckKffMatrix(self, kff)
            self.AFLAG = 1;
            if log10(condest(kff)) > 12
                self.AFLAG = 0;
            end 
            disp("Condition Number: " + round(condest(kff)));
            disp("Digits Lost: " + round(log10(condest(kff))));
        end
        
        %% Compute Displacements Reactions
        function [DEFL, REACT] = ComputeDisplacementsReactions(self, kff, ff, feff, kfn, ksf, ksn, nnodes)
            %disp deltas
            un = self.fix_t(self.dispDOFS);
            %free deltas
            uf = kff\((ff-feff) - kfn*un); % freexdisp * disp = free
            %reactions
            ps = ksf*uf + ksn*un;
            %empty matrix of zeros
            DEFL = zeros(6,nnodes);
            DEFL(self.freeDOFS) = uf;
            DEFL(self.dispDOFS) = un;
            DEFL = transpose(DEFL);
            REACT = zeros(6,nnodes);
            REACT(self.fixedDOFS) = ps;
            REACT = transpose(REACT); %might need to add reactions at disp dofs?
            disp("Displacements at free DOFS:")
            disp(uf)
        end
        
        %% Recover Element Forces
        function ELE_FOR = RecoverElementForces(self, DEFL)
            ELE_FOR = zeros(length(self.element_vector), 12);
            for i = 1:length(self.element_vector)
                ele_dofs = GetDOFS(self.element_vector(i));
                f_local = ComputeForces(self.element_vector(i), DEFL(ele_dofs));
                ELE_FOR(i,:) = transpose(f_local);
            end
        end
        
        %% Compute Error in Results
        function ComputeError(self, kff, kfn, ff, feff)
            transposedDEFL = transpose(self.DEFL);
            ff_computed = kff*transposedDEFL(self.freeDOFS);
            error = (ff - feff) - ff_computed - kfn*transposedDEFL(self.dispDOFS);
            disp("Error in Applied Load Vector");
            disp(error)
        end
    end
end
