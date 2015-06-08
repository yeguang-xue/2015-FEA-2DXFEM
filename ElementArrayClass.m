classdef ElementArrayClass < handle
    properties
        num_elements             % number of elements
        data =  cell(1)          % elements data
    end
    
    methods
        function obj=ElementArrayClass(n)
            if nargin ~= 0
                obj.num_elements = n;
                obj.data = cell(n,1);
            end
        end
        
        function obj=CreateElement(obj, id, element)
            obj.data{id} = element;
        end
        
        function element=getElementByID(obj, id)
            element = obj.data{id};
        end
        
        function plot(obj)
            for id = 1:obj.num_elements
                obj.data{id}.plot();
            end
        end
        
        function AssignMaterial(obj,material)
            for id = 1:obj.num_elements
                obj.data(id).AssignMaterial(material);
            end
        end
        
        function M = build_mass_matrix(obj,num_dof)
            M = zeros(num_dof);
            for id = 1:obj.num_elements
                element = obj.data(id);
                M = element.build_mass_matrix();
            end
        end
        
        function K = build_stiff_matrix(obj,num_dof)
            K = zeros(num_dof);
            for id = 1:obj.num_elements
                element = obj.getElementByID(id);
                conn = element.conn;
                
                Ke = element.build_stiff_matrix();
                doflist = obj.nodes.dofmapping(conn);
                K(doflist,doflist) =  K(doflist,doflist) + Ke;
            end
        end
        
        function CalcResults(obj)
            for id = 1:obj.num_elements
                element = obj.getElementByID(id);
                element.CalcResults();
            end
        end
        
        function outputVtk(obj,fid,outputtype)
            if strcmp(outputtype,'CELLS')
                fprintf(fid,'CELLS %d %d\n',obj.num_elements,obj.num_elements*5 );
                for id = 1:obj.num_elements
                    element = obj.getElementByID(id);
                    element.outputVtk(fid,'CELLS')
                end
            end
            if strcmp(outputtype,'CELL_TYPES')
                fprintf(fid,'CELL_TYPES %d\n',obj.num_elements);
                for id = 1:obj.num_elements
                    element = obj.getElementByID(id);
                    element.outputVtk(fid,'CELL_TYPES')
                end
            end
        end        
        
    end
end