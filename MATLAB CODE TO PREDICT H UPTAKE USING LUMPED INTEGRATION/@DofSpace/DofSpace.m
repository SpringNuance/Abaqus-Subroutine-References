classdef DofSpace < handle
    %DOFSPACE Management of degrees of freedom
    %   Class which keeps track of added degrees of freedom, giving index
	%   of dof for a dof name, and location for a combination of index and
	%   node.
    
    properties
        mesh		%pointer to mesh object
        DofTypes	%String vector containing names of each dof type
        DofNumbering   %node, dofType -> dofnumber
        NDofs		%total amout of degrees of freedom
    end
    
    methods
        function obj = DofSpace(mesh)
            %DOFSPACE Construct an instance of this class
            
            obj.mesh = mesh;
            obj.DofTypes = {};
            obj.NDofs = 0;
            obj.DofNumbering = sparse(length(obj.mesh.Nodes), length(obj.DofTypes));
        end
        
        function dofIndex = addDofType(obj, dofnames)
			% registers the name of the degree of freedom and returns its
			% index. If it already exists doesn't register it anew and
			% solely returns the index

            dofIndex = zeros(length(dofnames),1);
            dofTypeCount = length(obj.DofTypes);
            
            for i=1:length(dofnames)
               alreadyexist = false;
               for j=1:dofTypeCount
                  if (dofnames{i}==obj.DofTypes{j})
                      alreadyexist = true;
                      dofIndex(i) = j;
                  end
               end
               if (alreadyexist == false)
                   dofTypeCount = dofTypeCount+1;
                   obj.DofTypes{dofTypeCount} = dofnames{i};
                   dofIndex(i) = dofTypeCount;
                   obj.DofNumbering(:,length(obj.DofTypes)) = 0;
               end
			end
        end
        

        function addDofs(obj, dofIndices, nodeIndex)
			% Adds degrees of freedom for type "dofIndices" to the nodes
			% "nodeIndex"
            for i=1:length(dofIndices)
                for k=1:length(nodeIndex)
                    if (nodeIndex(k)<=size(obj.DofNumbering, 1))
                        curNum = obj.DofNumbering(nodeIndex(k), dofIndices(i));
                    else
                        obj.DofNumbering(nodeIndex(k), :) = 0;
                        curNum = 0;
                    end
                    if (curNum == 0) %dof does not yet exist if 0, add dof
                        obj.NDofs = obj.NDofs+1;
                        obj.DofNumbering(nodeIndex(k), dofIndices(i)) = obj.NDofs;
                    end
                end
            end
        end
        
        function DofTypeIndex = getDofType(obj, dofnames)
			% returns the dof type index for pre-existing degrees of
			% freedom
            DofTypeIndex = zeros(length(dofnames),1);
            for i=1:length(dofnames)
               for j=1:length(obj.DofTypes)
                  if (dofnames{i}==obj.DofTypes{j})
                      DofTypeIndex(i) = j;
                  end
               end
            end
        end
        
        function DofIndices = getDofIndices(obj, dofType, NodeIndices)
           % gets the indices for a combination of degree of freedom
		   % "doftype" and nodes "NodeIndices"
           DofIndices = obj.DofNumbering(NodeIndices, dofType);
           if (length(find(DofIndices==0)) ~= 0)
              disp("error here (dofspace)") 
           end
		end

    end
end

