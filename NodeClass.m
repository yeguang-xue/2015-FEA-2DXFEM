classdef NodeClass < handle
    properties
        X           % Reference coordinate
        u           % Displacements
        v           % Velocity
        q           % Node enrichment
        f           % Nodel force
        enrich = 0; % enrich = 0, no enrich; enrich = id, enrich id
        level       % Node level;
    end
    
    methods
        function obj=NodeClass(X)
            if nargin ~= 0
                obj.X = X;
            end
        end
        
        function Assign_coords(obj,X)
            obj.X = X;
        end
        
        function Assign_enrich(obj,enrich)
            obj.enrich = enrich;
        end
        
        function Assign_level(obj,level)
            obj.level = level;
        end
        
        function AssignResults(obj,d)
            if obj.enrich == 0
                obj.u = d;
            else
                obj.u = d([1,2]);
                obj.q = d([3,4]);
            end
        end
        
        function plot(obj)
            plot( obj.X(1), obj.X(2), 'ko')
        end
        
        
        function outputVtk(obj,fid,outputtype)
            
            if strcmp(outputtype,'POINTS_UNDEFORMED')
                fprintf(fid,'%f %f %f\n',[obj.X(1),obj.X(2),0]);
            end
            
            if strcmp(outputtype,'POINTS_DEFORMED')
                fprintf(fid,'%f %f %f\n',[obj.X(1)+obj.u(1),obj.X(2)+obj.u(2),0]);
            end
            
            if strcmp(outputtype,'U_DATA')
                fprintf(fid,'%f %f\n',[obj.u(1),obj.u(2)]);
            end
        end
    end
end