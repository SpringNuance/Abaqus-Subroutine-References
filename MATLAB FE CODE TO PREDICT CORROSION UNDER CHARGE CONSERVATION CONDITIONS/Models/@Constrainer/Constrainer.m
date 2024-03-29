classdef Constrainer < BaseModel
    %Constrains specified degrees of freedom, input example
% 		physics_in{4}.type = "Constrainer";   %name of this model
%		physics_in{4}.Egroup = "M_Bottom";		%name of element group to constrain
%		physics_in{4}.dofs = {"dy"};			%name of degree of freedom to constrain
%		physics_in{4}.conVal = [0];				%value it should be constrained to
    
    properties
        myName			%name of this model
        mesh			%pointer to mesh object
        myGroup			%constrinaed element group name
        myGroupIndex	%element group number
        dofSpace		%pointer towards the dofspace
        dofTypeIndices	%index of constrained dofs
        
        conVal			%value to constrain degrees of freedom to
    end
    
    methods
        function obj = Constrainer(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Constrainer";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getNodeGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType(inputs.dofs);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForNodeGroup(obj.myGroupIndex));
            
            obj.conVal = inputs.conVal;
        end
        
        function getKf(obj, physics)
			% saves constraints to constraint object contained within the
			% physics. This does not directly apply the constraints, they
			% are merely collected, after which the physics object does the
			% actual contraining 

            allNodes = obj.mesh.GetAllNodesForNodeGroup(obj.myGroupIndex);
            for i=1:length(obj.dofTypeIndices)
                newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i), allNodes);
                
                physics.condofs = [physics.condofs; newcons];
                physics.convals = [physics.convals; newcons*0+obj.conVal(i)];
            end
        end
    end
end

