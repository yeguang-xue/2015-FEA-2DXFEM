classdef CPS4XElementClass < ElementClass
    properties
        material_1               % Material 1
        material_2               % Material 2
        gauss                    % Gauss position
        weight                   % Weight
        quad                     % quadrature points
        vertex_x = zeros(4,1);
        vertex_y = zeros(4,1);
        vertex_level = zeros(4,1);   % level at nodes
        strain                   %
    end
    
    methods (Static)
        function Nlist = shape_parent(X)
            xi = X(1);
            eta = X(2);
            
            N1 = 1/4*(1-xi)*(1-eta);
            N2 = 1/4*(1+xi)*(1-eta);
            N3 = 1/4*(1+xi)*(1+eta);
            N4 = 1/4*(1-xi)*(1+eta);
            
            Nlist = [ N1, N2, N3, N4];
        end
        
        function Blist = derative_parent(X)
            xi = X(1);
            eta = X(2);
            
            B11 = - 1/4*(1-eta);
            B12 = - 1/4*(1-xi);
            B21 =   1/4*(1-eta);
            B22 = - 1/4*(1+xi);
            B31 =   1/4*(1+eta);
            B32 =   1/4*(1+xi);
            B41 = - 1/4*(1+eta);
            B42 =   1/4*(1-xi);

            Blist = [ B11, B21, B31, B41;
                B12, B22, B32, B42];
        end
    end
    
    methods
        
        function obj=CPS4XElementClass(nodes,conn,material_1,material_2,level)
            if nargin ~= 5
                nodes = NodeArrayClass;
                conn = [];
                material_1 = [];
                material_2 = [];
                level = [];
            end
            
            obj.type = 'CPS4X';
            obj.nodes = nodes;
            obj.conn = conn;
            obj.material_1 = material_1;
            obj.material_2 = material_2;
            obj.vertex_level = level;
            
            % Assign enrich information
            for kk = 1:4
                node = obj.nodes.getNodeByID(conn(kk));
                node.Assign_enrich(1);
                node.Assign_level(level(kk));
            end
                
            % Calculate vertex_x, vertex_y
            for kk = 1:4
                node = obj.nodes.getNodeByID(conn(kk));
                obj.vertex_x(kk) = node.X(1);
                obj.vertex_y(kk) = node.X(2);
            end
            
            xn = obj.vertex_x;
            yn = obj.vertex_y;
                        
            % Get gauss points and weights
            [obj.gauss,obj.weight] = getGaussPoints(6);
            
            % Calculate quadrature points postions
            for kk = 1:size(obj.gauss,1)
                xi  = obj.gauss(kk,1);
                eta = obj.gauss(kk,2);
                
                Nlist = obj.shape_parent([xi,eta]);
                
                obj.quad(kk,1) = Nlist*xn;
                obj.quad(kk,2) = Nlist*yn;
            end
        end
        
        function plot(obj)
            plot( obj.vertex_x([1:end,1]), obj.vertex_y([1:end,1]), '-')
%             for kk = 1:size(obj.gauss,1)
%                 level = obj.local_ls(obj.gauss(kk,:));
%                 if level>0
%                     plot( obj.quad(kk,1), obj.quad(kk,2), 'b*')
%                 elseif level<0
%                     plot( obj.quad(kk,1), obj.quad(kk,2), 'r*')
%                 else
%                     plot( obj.quad(kk,1), obj.quad(kk,2), 'k*')
%                 end
%             end
        end
        
        function Je = jacob_fun(obj,X)
            xn = obj.vertex_x;
            yn = obj.vertex_y;
            Blist = obj.derative_parent(X);
            Je = Blist*[xn,yn];
        end
        
        function level = local_ls(obj,X)
            level = 0;
            Nlist = obj.shape_parent(X);
            for kk = 1:4
                level = level + Nlist(kk)*obj.vertex_level(kk);
            end
        end
        
        function PSI = enrich_fun(obj,X)
            PSI = 0;
            Nlist = obj.shape_parent(X);
            for kk = 1:4
                PSI = PSI + Nlist(kk)*obj.vertex_level(kk);
            end
            PSI =abs(PSI);
        end
        
        function eX_Blist = eX_derative_parent(obj,X)
            
            Nlist = obj.shape_parent(X);
            Blist = obj.derative_parent(X);
            PSI = obj.enrich_fun(X);                        
            level = obj.local_ls(X);
            
            PSI_xi = 0;
            PSI_eta = 0;
            
            if level >0
                for kk = 1:4
                    PSI_xi = PSI_xi + Blist(1,kk)*obj.vertex_level(kk);
                    PSI_eta = PSI_eta + Blist(1,kk)*obj.vertex_level(kk);
                end
            elseif level<0
                for kk = 1:4
                    PSI_xi = PSI_xi - Blist(1,kk)*obj.vertex_level(kk);
                    PSI_eta = PSI_eta - Blist(1,kk)*obj.vertex_level(kk);
                end
            end
            
            eX_Blist = zeros(2,4);
            
            for kk = 1:4 
                eX_Blist(1,kk) = Blist(1,kk)*PSI+PSI_xi*Nlist(kk);
                eX_Blist(2,kk) = Blist(2,kk)*PSI+PSI_eta*Nlist(kk);
            end
        end
        
        function Be = derative_matrix(obj,X)
            Blist = obj.derative_parent(X);
            B = obj.jacob_fun(X)\Blist;
            B = B';
            
            Be = [B(1,1),      0, B(2,1),      0, B(3,1),      0, B(4,1),      0;
                       0, B(1,2),      0, B(2,2),      0, B(3,2),      0, B(4,2);
                  B(1,2), B(1,1), B(2,2), B(2,1), B(3,2), B(3,1), B(4,2), B(4,1)];
        end
        
        function eX_Be = eX_derative_matrix(obj,X)
            eX_Blist = obj.eX_derative_parent(X);
            eX_B = obj.jacob_fun(X)\eX_Blist;
            eX_B = eX_B';
            
            eX_Be = [eX_B(1,1),      0, eX_B(2,1),      0, eX_B(3,1),      0, eX_B(4,1),      0;
                       0, eX_B(1,2),      0, eX_B(2,2),      0, eX_B(3,2),      0, eX_B(4,2);
                  eX_B(1,2), eX_B(1,1), eX_B(2,2), eX_B(2,1), eX_B(3,2), eX_B(3,1), eX_B(4,2), eX_B(4,1)];
        end
        
        function Me = build_mass_matrix(obj)
            rho = obj.material(1);
            Me = rho;
        end
        
        function Ke = build_stiff_matrix(obj)
            youngs_1 = obj.material_1(2);
            nu_1 = obj.material_1(3);
            youngs_2 = obj.material_2(2);
            nu_2 = obj.material_2(3);
            
            De_1 = [1, nu_1, 0; nu_1, 1, 0; 0, 0, (1-nu_1)/2];
            De_1 = De_1*youngs_1/(1-nu_1^2);
            De_2 = [1, nu_2, 0; nu_2, 1, 0; 0, 0, (1-nu_2)/2];
            De_2 = De_2*youngs_2/(1-nu_2^2);
            
            Kuu = zeros(8,8);
            Kuq = zeros(8,8);
            Kqq = zeros(8,8);
            
            for kk = 1:size(obj.gauss,1)
                level = obj.local_ls(obj.gauss(kk,:));
                Je = obj.jacob_fun(obj.gauss(kk,:));
                Be = obj.derative_matrix(obj.gauss(kk,:));
                eX_Be = obj.eX_derative_matrix(obj.gauss(kk,:));
                Je_det = det(Je);
                if level>0
                    Kuu = Kuu + obj.weight(kk)*Be'*De_1*Be*Je_det;
                    Kuq = Kuq + obj.weight(kk)*Be'*De_1*eX_Be*Je_det;
                    Kqq = Kqq + obj.weight(kk)*eX_Be'*De_1*eX_Be*Je_det;
                else
                    Kuu = Kuu + obj.weight(kk)*Be'*De_2*Be*Je_det;
                    Kuq = Kuq + obj.weight(kk)*Be'*De_2*eX_Be*Je_det;
                    Kqq = Kqq + obj.weight(kk)*eX_Be'*De_2*eX_Be*Je_det;
                end
            end
            
            Ke = [Kuu, Kuq;Kuq', Kqq];
        end
        
        function CalcResults(obj)
            obj.strain = 0;
        end
        
        function outputVtk(obj,fid,outputtype)
            if strcmp(outputtype,'CELLS')
                fprintf(fid,'%d %d %d %d %d\n',[4, obj.conn-1]);
            end
            if strcmp(outputtype,'CELL_TYPES')
                fprintf(fid,'%d\n', 9);
            end
        end
    end
end