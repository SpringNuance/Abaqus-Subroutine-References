classdef BaseModel < handle
    %Base model from which all other models are derived (solely functions
	%to prevent having to declare the same function multiple times)
    
    properties
        
    end
    
    methods
        function obj = BaseModel()

        end
        
        function Commit(obj, physics, commit_type)
            
        end
        
        function Irr = Irreversibles(obj, physics)
            Irr = false;
        end
        
        function [hasInfo, provided] = Provide_Info(obj, physics, var, elems, loc)
           hasInfo = false;
           provided = [];
        end
        
        function getKf(obj, physics)

        end
    end
end

