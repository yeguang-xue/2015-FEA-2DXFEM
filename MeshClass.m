classdef MeshClass < handle
    properties
        nodes =  NodeArrayClass         %
        elements =  ElementArrayClass   %
        num_dof                         %
    end
    
    methods
        function obj = MeshClass(nodes, elements)
            obj.nodes = nodes;
            obj.elements = elements;
            obj.num_dof = 0;
            enrichid = 0;
            
            for id = 1:obj.nodes.num_nodes
                node = obj.nodes.getNodeByID(id);
                if node.enrich == 0
                    obj.num_dof = obj.num_dof +2;
                else
                    enrichid = enrichid +1;
                    node.Assign_enrich(enrichid);
                    obj.num_dof = obj.num_dof +4;
                end
            end
        end
        
        function plot(obj)
            obj.elements.plot()
            obj.nodes.plot()
        end
        
        function M = build_mass_matrix(obj)
            M = obj.elements.build_mass_matrix(obj.num_dof);
        end
        
        function K = build_stiff_matrix(obj)
            K = zeros(obj.num_dof);
            for id = 1:obj.elements.num_elements
                element = obj.elements.getElementByID(id);
                conn = element.conn;
                Ke = element.build_stiff_matrix();
                doflist = element.dofmapping(conn);
                K(doflist,doflist) =  K(doflist,doflist) + Ke;
            end
        end
        
        function f = build_force_matrix(obj,forceBC)
            f = zeros(obj.num_dof,1);
            for fid = 1:size(forceBC,1)
                nid = forceBC(bcid,1);
                dof = obj.nodes.dofmapping(nid);
                for dir = 1:2
                    fvalue = dispBC(bcid,dir+1);
                    f(dof(dir)) = fvalue;
                end
            end
        end
        
        function [K_modify, f_modify] = applydispBC(obj,K,f,dispBC)
            K_modify = K;
            f_modify = f;
            for bcid = 1:size(dispBC,1)
                nid = dispBC(bcid,1);
                dof = obj.nodes.dofmapping(nid);
                for dir = 1:2
                    uvalue = dispBC(bcid,dir+1);
                    if ~isnan(uvalue)
                        K_modify(dof(dir),:) = 0;
                        K_modify(dof(dir),dof(dir)) = 1;
                        f_modify(dof(dir)) = uvalue;
                    end
                end
            end
        end
        
        function AssignResults(obj,d)
            obj.nodes.AssignResults(d);
            obj.elements.CalcResults();
        end
        
        function outputVtk(obj,filename,outcontrol)
            
            if nargin <= 2
                outcontrol = 'UNDEFORMED';
            end
            
            fid = fopen(filename,'w');
            fprintf(fid,'# vtk DataFile Version 2.0\n');
            fprintf(fid,'Results from FEA\n');
            fprintf(fid,'ASCII\n');
            
            fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
            
            if strcmp(outcontrol,'DEFORMED')
                obj.nodes.outputVtk(fid,'POINTS_DEFORMED');
            else
                obj.nodes.outputVtk(fid,'POINTS_UNDEFORMED');
            end
            
            obj.elements.outputVtk(fid,'CELLS');
            obj.elements.outputVtk(fid,'CELL_TYPES');      
            
            obj.nodes.outputVtk(fid,'POINT_DATA');
            
            fclose(fid);
        end
    end
end