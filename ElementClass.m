%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yeguang Xue (Northwestern Univ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ElementClass < handle
    properties
        type                     % Element type
        nodes                    % Nodes list
        conn                     % Connect list
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
        function AssignNodes(obj,nodes)
            obj.nodes = nodes;
        end
        
        function AssignConn(obj,conn)
            obj.conn = conn;
        end
        
        function AssignMaterial(obj,material)
            obj.material = material;
        end
        
        function doflist = dofmapping(obj,id)
            doflist = [];
            udof = [];
            qdof = [];
            for kk = 1:length(id)
                node = obj.nodes.getNodeByID(id(kk));
                tmpudof = [id(kk)*2-1, id(kk)*2];
                udof = [udof, tmpudof];
                if node.enrich ~= 0
                    tmpqdof = [obj.nodes.num_nodes*2+node.enrich*2-1, obj.nodes.num_nodes*2+node.enrich*2];
                    qdof = [qdof, tmpqdof];
                end
            end

            if strcmp(obj.type,'CPS4')
                doflist = udof;
            elseif strcmp(obj.type,'CPS4X')
                doflist = [udof,qdof];
            end 
        end
        
    end
    
    methods (Abstract)
        plot(obj)
        % plotresult(obj)
        %         N = shape_fun(obj,X)
        %         Be = derative_matrix(obj,X)
        Me = build_mass_matrix(obj)
        Ke = build_stiff_matrix(obj)
    end
end