classdef Constrainer < BaseModel
    %CONSTRAINER Physics model that applies constraints to the governing
	%equations, based on a provided element group. Required inputs:
    %	physics_in{3}.type = "Constrainer";		Name of this class
    %	physics_in{3}.Egroup = "M_Bottom";		Element group which is being constrained
    %	physics_in{3}.dofs = {"dx"};			Name of the Degree of Freedom being constrained
    %	physics_in{3}.conVal = [0];				Value of the constraint being applied
    
    properties
        myName			%Name of this model used for identification purposes
        mesh			%Pointer to the mesh object
        myGroup			%String to indicate the name of the element group being constrained
        myGroupIndex	%index of the constrained group
        dofSpace		%Pointer to the degree of freedom object
        dofTypeIndices	%Indices referring to the constrained dofs
        conVal			%Values that are applied as constraint
    end
    
    methods
        function obj = Constrainer(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Constrainer";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType(inputs.dofs);
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
            
            obj.conVal = inputs.conVal;
        end
        
        function getKf(obj, physics)
			% saves constraints to constraint object contained within the
			% physics. This does not directly apply the constraints, they
			% are merely collected, after which the physics object does the
			% actual contraining 

            allNodes = obj.mesh.GetAllNodesForGroup(obj.myGroupIndex);
            for i=1:length(obj.dofTypeIndices)
                newcons = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i), allNodes);
                
                physics.condofs = [physics.condofs; newcons];
                physics.convals = [physics.convals; newcons*0+obj.conVal(i)];
            end
        end
    end
end

