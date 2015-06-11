%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yeguang Xue (Northwestern Univ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef CPS4ElementClass < ElementClass
    properties
        material                 % Material
        gauss                    % Gauss position
        weight                   % Weight
        quad                     % quadrature points
        vertex_x = zeros(4,1);
        vertex_y = zeros(4,1);
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
        
        function obj=CPS4ElementClass(nodes,conn,material)
            if nargin ~= 3
                nodes = NodeArrayClass;
                conn = [];
                material = [];
            end
            
            obj.type = 'CPS4';
            obj.nodes = nodes;
            obj.conn = conn;
            obj.material = material;
            
            % Calculate vertex_x, vertex_y
            for kk = 1:4
                node = obj.nodes.getNodeByID(conn(kk));
                obj.vertex_x(kk) = node.X(1);
                obj.vertex_y(kk) = node.X(2);
            end
            
            xn = obj.vertex_x;
            yn = obj.vertex_y;
                        
            % Get gauss points and weights
            [obj.gauss,obj.weight] = getGaussPoints(2);
            
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
            plot( obj.vertex_x([1:end,1]), obj.vertex_y([1:end,1]), 'k-')
            plot( obj.quad(1:end,1), obj.quad(1:end,2), '*')
            if obj.material(2)>2
                fill( obj.vertex_x([1:end,1]), obj.vertex_y([1:end,1]), 'r')
            else
                fill( obj.vertex_x([1:end,1]), obj.vertex_y([1:end,1]), 'b')
            end
        end
        
        function Je = jacob_fun(obj,X)
            xn = obj.vertex_x;
            yn = obj.vertex_y;
            Blist = obj.derative_parent(X);
            Je = Blist*[xn,yn];
        end
        
        function Be = derative_matrix(obj,X)
            Blist = obj.derative_parent(X);
            B = obj.jacob_fun(X)\Blist;
            B = B';
            
            Be = [B(1,1),      0, B(2,1),      0, B(3,1),      0, B(4,1),      0;
                       0, B(1,2),      0, B(2,2),      0, B(3,2),      0, B(4,2);
                  B(1,2), B(1,1), B(2,2), B(2,1), B(3,2), B(3,1), B(4,2), B(4,1)];
        end
        
        
        function Me = build_mass_matrix(obj)
            rho = obj.material(1);
            Me = rho;
        end
        
        function Ke = build_stiff_matrix(obj)
            youngs = obj.material(2);
            nu = obj.material(3);
            
            De = [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
            De = De*youngs/(1-nu^2);
            
            Ke = zeros(8,8);
            
            for kk = 1:size(obj.gauss,1)
                Je = obj.jacob_fun(obj.gauss(kk,:));
                Be = obj.derative_matrix(obj.gauss(kk,:));
                Je_det = det(Je);
                Ke = Ke + obj.weight(kk)*Be'*De*Be*Je_det;
            end
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