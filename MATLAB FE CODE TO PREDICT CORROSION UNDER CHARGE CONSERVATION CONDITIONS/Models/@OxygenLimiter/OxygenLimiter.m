classdef OxygenLimiter < BaseModel
    %Constrains specified degrees of freedom up to a set time, input example
	%	physics_in{3}.type = "OxygenLimiter";
	%	physics_in{3}.Egroup = "E_Top";
	%	physics_in{3}.dofs = {"O2"}; %Name of degree of freedom
	%	physics_in{3}.conVal = [initO2];	%value of oxygen constraint
	%	physics_in{3}.tmax = 48*3600;	%time after which constraint is removed
    
    properties
        myName			%name of this model
        mesh			%pointer to mesh object
        myGroup			%constrinaed element group name
        myGroupIndex	%element group number
        dofSpace		%pointer towards the dofspace
        dofTypeIndices	%index of constrained dofs
        
        conVal			%value to constrain degrees of freedom to
		tmax			%Time untill constraint is removed
    end
    
    methods
        function obj = OxygenLimiter(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "OxygenLimiter";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getNodeGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType(inputs.dofs);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForNodeGroup(obj.myGroupIndex));
            
            obj.conVal = inputs.conVal;
			obj.tmax = inputs.tmax;
        end
        
        function getKf(obj, physics)
			time = physics.time;
			if (time<obj.tmax)
            	allNodes = obj.mesh.GetAllNodesForNodeGroup(obj.myGroupIndex);
            	
				for i=1:length(obj.dofTypeIndices)
                	newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i), allNodes);
                	
                	physics.condofs = [physics.condofs; newcons];
                	physics.convals = [physics.convals; newcons*0+obj.conVal(i)];
				end
			end
        end
    end
end

