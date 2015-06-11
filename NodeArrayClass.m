%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yeguang Xue (Northwestern Univ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef NodeArrayClass < handle
    properties
        num_nodes          % number of nodes
        data = NodeClass   % nodes data
    end
    
    methods
        function obj=NodeArrayClass(n)
            if nargin ~= 0
                obj.num_nodes = n;
                obj.data(n) = NodeClass;
            end
        end
        
        function CreateNode(obj, id, node)
            obj.data(id) = node;
        end
        
        function node=getNodeByID(obj, id)
            node = obj.data(id);
        end
        
        function dof = dofmapping(obj,id)
            node = obj.getNodeByID(id);
            if node.enrich == 0
                tmpdof = [id*2-1, id*2];
            else
                udof = [id*2-1, id*2];
                qdof = [obj.num_nodes*2+node.enrich*2-1, obj.num_nodes*2+node.enrich*2];
                tmpdof = [udof, qdof];
            end
            dof = tmpdof;
        end
        
        function AssignResults(obj,d)
            for id = 1:obj.num_nodes
                node = obj.getNodeByID(id);
                dof = obj.dofmapping(id);
                nodevalue = d(dof);
                node.AssignResults(nodevalue);
            end
        end
        
        function plot(obj)
            for id = 1:obj.num_nodes
                node = obj.getNodeByID(id);
                node.plot();
            end
        end
        
        function outputVtk(obj,fid,outputtype)
            % DEFORMED NODES POSITIONS
            if strcmp(outputtype,'POINTS_DEFORMED')
                fprintf(fid,'POINTS %d float\n',obj.num_nodes);
                for id = 1:obj.num_nodes
                    node = obj.getNodeByID(id);
                    node.outputVtk(fid,'POINTS_DEFORMED');
                end
            end
            % UNDEFORMED NODES POSITIONS            
            if strcmp(outputtype,'POINTS_UNDEFORMED')
                fprintf(fid,'POINTS %d float\n',obj.num_nodes);
                for id = 1:obj.num_nodes
                    node = obj.getNodeByID(id);
                    node.outputVtk(fid,'POINTS_UNDEFORMED');
                end
            end
            % NODAL DATA
            if strcmp(outputtype,'POINT_DATA')
                fprintf(fid,'POINT_DATA %d\n',obj.num_nodes);
                
                fprintf(fid,'FIELD FieldData 1\n');
                fprintf(fid,'U 2 %d float\n',obj.num_nodes);
                
                for id = 1:obj.num_nodes
                    node = obj.getNodeByID(id);
                    node.outputVtk(fid,'U_DATA');
                end
            end
        end
    end
end