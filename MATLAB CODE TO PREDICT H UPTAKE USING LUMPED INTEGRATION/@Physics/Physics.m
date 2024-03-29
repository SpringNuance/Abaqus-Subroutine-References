classdef Physics < handle
    %PHYSICS Class that handles the definitions and assembly of stiffness
	%matrices, constraints, and statevectors
    
    properties
        mesh	%pointer to the mesh object
        models	%array containing all the included physical models
        dofSpace %pointer to the degree of freedom object
        
        K	%Tangential stiffness matrix
        fint	%Internal force vector
        StateVec	%State vector which is iteratively updated to obtain results at t+dt
        StateVec_Old	%Converged state vector at time t
        time	%current time within the simulation
        dt	%Time increment size
        
        condofs	%constrained degrees of freedom
        convals %values applied as constraints to the degrees of freedom
        convals_corr	%increment required to maintain constraint
        conMat		%matrix to transfer from unconstrained system to solely the constrained dofs
        unconMat	%matrix to transfer from unconstrained system to constrained
    end
    
    methods
        PlotNodal(obj, dofName, dispscale, plotloc) %exterior defined, plots nodal quantities
        PlotIP(obj, varName, plotloc) %exterior defined, plots integration point quantities
        
        function obj = Physics(mesh, inputs)
			%initialization
            obj.mesh = mesh;
            obj.dofSpace = DofSpace(obj.mesh);
            
            for i=1:length(inputs)
                f = str2func(inputs{i}.type);
                obj.models{i} = f(mesh, obj, inputs{i});
            end
            
            obj.dt = 0;
            
            dofcount = obj.dofSpace.NDofs;
            obj.StateVec = zeros(dofcount, 1);
            obj.StateVec_Old = obj.StateVec;
            
            obj.K = sparse(dofcount, dofcount);
            obj.fint = zeros(dofcount,1);
        end
        
        function Assemble(obj)
			%Assemble stiffness matrix and internal force vector
            dofcount = obj.dofSpace.NDofs;

			obj.condofs = [];
			obj.convals = [];

            nonz = round(nnz(obj.K)*1.2);
            obj.K = spalloc(dofcount, dofcount, nonz);
            obj.fint = zeros(dofcount, 1);

            disp("    Assembling:")
            for m=1:length(obj.models)
                obj.models{m}.getKf(obj);
            end
        end
       
        function Commit(obj, commit_type)
			% commit time dependent state on progressing to the next time
			% increment
            for m=1:length(obj.models)
                obj.models{m}.Commit(obj, commit_type);
            end
            
            if (commit_type == "Timedep")
                obj.StateVec_Old = obj.StateVec;
            end
        end
        
        function anyIrr = Irreversibles(obj)
			% Checks if any irreversable processes can occur at the end of
			% the time step 
            anyIrr = false;
            for m=1:length(obj.models)
                anyIrr = anyIrr + obj.models{m}.Irreversibles(obj);
            end
        end
            
        function Constrain(obj)
			%constrain the tangential stiffness matrix and force vector,
			%based on defined constraints

            obj.convals_corr = obj.convals - obj.StateVec(obj.condofs);
            basemat = speye(size(obj.K));
            obj.unconMat = basemat;
            obj.unconMat(:, obj.condofs) = [];
            obj.conMat = basemat(:, obj.condofs);

            obj.fint = obj.unconMat'*obj.fint + obj.unconMat'*obj.K*obj.conMat*obj.convals_corr;
            obj.K    = obj.unconMat'*obj.K*obj.unconMat;
        end
        
        function Update(obj, dx)
			% adds increment dx to the state vector
            obj.StateVec = obj.StateVec + obj.unconMat*dx + obj.conMat*obj.convals_corr;
        end
        
        function info = Request_Info(obj, var, elems, loc)
			%inter-model communicator able to request information based on
			%variable name var
            info = false;
            for m=1:length(obj.models)
                [hasInfo, provided] = obj.models{m}.Provide_Info(obj, var, elems, loc);
                if (hasInfo)
                    info = provided;
                end
            end   
        end

    end
end

