classdef Electrolyte < BaseModel
    %Model which uses the nernst-planck equations to describe the diffusion
	%of ions within a non-flowing electrolyte
    %   inputs example: 
	% 		physics_in{1}.type = "Electrolyte";   %name of this model
	% 	    physics_in{1}.Egroup = "Electrolyte"; %domain name on which the model works
	% 	    physics_in{1}.D = [9.3; 5.3; 1.3; 2; 1.4; 1]*1e-9;  %H OH Na CL Fe FeOH, diffusion coefficient of ion species
	% 	    physics_in{1}.z = [1; -1; 1; -1; 2; 1]; %ion charges
	% 	    physics_in{1}.pH0 = 5;  %initial pH
	%		physics_in{1}.O2 = initO2; %initial O2 concentration
	% 	    physics_in{1}.NaCl = 0.6e3;  %initial concentration of Na and CL
	% 	    physics_in{1}.Lumped = [true; true]; %water, metal; Should lumped integration be used for these volume reactions
	% 	    physics_in{1}.k = [1e6; 1e-1; 1e-3; 1e-3]; %water, Fe, Fe', FeOH; reaction rates for volume reactions
    
    properties
        mesh			%pointer towards the mesh
        myName			%name of this model
        myGroup			%element group on which it operates
        myGroupIndex	%element group index
        dofSpace		%pointer towards the dofspace
        dofTypeIndices	%indices of added degrees of freedom

		DoOnce			%flag to initialise some variables first time this runs

        D		%diffusion coefficient vector
		z		%charge vector
		k		%reaction rate vector
		Lumped	%flag vector for whether to use lumped integration
		pH0		%initial pH
		O2_0    %initial O2
		NaCl	%initial filler ion concentrations

		F_const = 96485.3329;		%faraday constant
		R_const = 8.31446261815324;	%gas constant
		T_const = 293.15;			%environment temperature

		n_species					%total number of ion species simulated
    end
    
    methods
        function obj = Electrolyte(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "Electrolyte";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.myGroup = inputs.Egroup;
            obj.myGroupIndex = obj.mesh.getGroupIndex(obj.myGroup);
            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs
            obj.dofTypeIndices = obj.dofSpace.addDofType({"Epot","H", "OH","Na", "Cl","Fe","FeOH", "O2"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));

			obj.DoOnce = true;
            
            %% get parameters
            obj.D = inputs.D;
			obj.pH0 = inputs.pH0;
			obj.NaCl = inputs.NaCl;
			obj.O2_0 = inputs.O2;
			obj.z = inputs.z;
			obj.k=inputs.k;
			obj.Lumped=inputs.Lumped;

			obj.n_species = length(obj.D);
        end
        
        function getKf(obj, physics)
			% Force vector and tangential matrix assembly procedure
            fprintf("        Electrolyte get Matrix:")
            t = tic;
            
            dt = physics.dt;

			if (obj.DoOnce) %set-up required to be erformed once
				%set initial values
				dofsE = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
				for i=1:obj.n_species
					dofsC(:,i) = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i+1), obj.mesh.GetAllNodesForGroup(obj.myGroupIndex));
				end
				initH = 1000*10^(-obj.pH0);
				initOH = 1000*10^(-14+obj.pH0);
				initCl = obj.NaCl;
				initNa = initCl-initH+initOH;
				initO2 = obj.O2_0;

				initState = [initH, initOH, initNa, initCl, 0, 0, initO2];

				physics.StateVec(dofsE) = 0;	physics.StateVec_Old(dofsE) = 0;
				for i=1:obj.n_species
					physics.StateVec(dofsC(:,i)) = initState(i);	
					physics.StateVec_Old(dofsC(:,i)) = initState(i);
				end
				obj.DoOnce = false;
			end
			clear dofsC dofsE

			%stiffness matrix
            dofmatX = [];
            dofmatY = [];
            kmat = [];
            fvec = [];
            dofvec = [];

			Svec = physics.StateVec;
            SvecOld = physics.StateVec_Old;
            parfor n_el=1:(size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)) %parfor
				% get nodes and element shaoe functions
                Elem_Nodes = obj.mesh.getNodes(obj.myGroupIndex, n_el);
                [N, G, w] = obj.mesh.getVals(obj.myGroupIndex, n_el);
				ipcoords = obj.mesh.getIPCoords(obj.myGroupIndex, n_el);
				w = w.*2.*pi.*ipcoords(1,:);

				%get nodal values
				dofsE = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
				dofsC = zeros(length(dofsE),obj.n_species);
				for i=1:obj.n_species
					dofsC(:,i) = obj.dofSpace.getDofIndices(obj.dofTypeIndices(i+1), Elem_Nodes);
				end

				C = Svec(dofsC);
				E = Svec(dofsE);

				C_OLD = SvecOld(dofsC);
				E_OLD = SvecOld(dofsE);

				%initialize arrays
				q_C    = zeros(length(dofsE), obj.n_species);
				dqC_dC = zeros(length(dofsE), length(dofsE), obj.n_species, obj.n_species);
				dqC_dE = zeros(length(dofsE), length(dofsE), obj.n_species);

				q_E = zeros(length(dofsE), 1);	
				dqE_dE = zeros(length(dofsE), length(dofsE));
				dqE_dC = zeros(length(dofsE), length(dofsE), obj.n_species);

				C_Lumped = zeros(length(dofsE), 1);
				for ip=1:length(w)

					%capacity
				%	for s=1:obj.n_species
				%		q_C(:,s) = q_C(:,s) + w(ip)/dt*N(ip,:)'*N(ip,:)*(C(:,s)-C_OLD(:,s));
				%		dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + w(ip)/dt*N(ip,:)'*N(ip,:);
				%	end
	
					%diffusion
					for s=1:obj.n_species
						Df = obj.D(s);

						q_C(:,s) = q_C(:,s) + w(ip)*squeeze(G(ip,:,:))*Df*squeeze(G(ip,:,:))'*C(:,s);
						dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + w(ip)*squeeze(G(ip,:,:))*Df*squeeze(G(ip,:,:))';
					end
	
					%electroneutrality
					for s=1:obj.n_species
						Df = obj.D(s);

						fct = obj.z(s)*obj.F_const/obj.R_const/obj.T_const*Df;
						
						q_C(:,s) =	q_C(:,s) + w(ip)*(N(ip,:)*C(:,s))*squeeze(G(ip,:,:))*(fct*squeeze(G(ip,:,:))'*E);		
						dqC_dC(:,:,s,s) = dqC_dC(:,:,s,s) + w(ip)*squeeze(G(ip,:,:))*(fct*squeeze(G(ip,:,:))'*E)*N(ip,:);
						dqC_dE(:,:,s) = dqC_dE(:,:,s) + w(ip)*(N(ip,:)*C(:,s))*squeeze(G(ip,:,:))*fct*squeeze(G(ip,:,:))';

						q_E = q_E + w(ip)*obj.z(s)*N(ip,:)'*N(ip,:)*C(:,s);
						dqE_dC(:,:,s) = dqE_dC(:,:,s) + w(ip)*obj.z(s)*N(ip,:)'*N(ip,:);
					end

					%% volume reactions

					% water
					if (obj.Lumped(1)==false)
						rwc = obj.k(1);
						rw = rwc*(1e-8-(N(ip,:)*C(:,1))*(N(ip,:)*C(:,2)));
						drw_dH  = rwc*(-(N(ip,:))*(N(ip,:)*C(:,2)));
						drw_dOH = rwc*(-(N(ip,:)*C(:,1))*(N(ip,:)));
						q_C(:,1) = q_C(:,1) - w(ip)*N(ip,:)'*rw;
						q_C(:,2) = q_C(:,2) - w(ip)*N(ip,:)'*rw;
	
						dqC_dC(:,:,1,1) = dqC_dC(:,:,1,1) - w(ip)*N(ip,:)'*drw_dH;
						dqC_dC(:,:,1,2) = dqC_dC(:,:,1,2) - w(ip)*N(ip,:)'*drw_dOH;
						dqC_dC(:,:,2,1) = dqC_dC(:,:,2,1) - w(ip)*N(ip,:)'*drw_dH;
						dqC_dC(:,:,2,2) = dqC_dC(:,:,2,2) - w(ip)*N(ip,:)'*drw_dOH;
					end

					%Fe
					if (obj.Lumped(2)==false)
						rcc = [obj.k(2), obj.k(3), obj.k(4)];
						rc = rcc(1)*(N(ip,:)*C(:,5))-rcc(2)*(N(ip,:)*C(:,6))*(N(ip,:)*C(:,1));
						drc_dH = -rcc(2)*(N(ip,:)*C(:,6))*N(ip,:);
						drc_dFE = rcc(1)*N(ip,:);
						drc_dFEOH = -rcc(2)*(N(ip,:)*C(:,1))*N(ip,:);
	
						rc2 = rcc(3)*(N(ip,:)*C(:,6));
						drc2_dFEOH = rcc(3)*N(ip,:);

						q_C(:,1) = q_C(:,1) - w(ip)*N(ip,:)'*(rc+rc2);
						q_C(:,5) = q_C(:,5) - w(ip)*N(ip,:)'*-rc;
						q_C(:,6) = q_C(:,6) - w(ip)*N(ip,:)'*(rc-rc2);
	
						dqC_dC(:,:,1,1) = dqC_dC(:,:,1,1) - w(ip)*N(ip,:)'*drc_dH;
						dqC_dC(:,:,1,5) = dqC_dC(:,:,1,5) - w(ip)*N(ip,:)'*drc_dFE;
						dqC_dC(:,:,1,6) = dqC_dC(:,:,1,6) - w(ip)*N(ip,:)'*(drc_dFEOH+drc2_dFEOH);
	
						dqC_dC(:,:,5,1) = dqC_dC(:,:,5,1) - w(ip)*N(ip,:)'*-drc_dH;
						dqC_dC(:,:,5,5) = dqC_dC(:,:,5,5) - w(ip)*N(ip,:)'*-drc_dFE;
						dqC_dC(:,:,5,6) = dqC_dC(:,:,5,6) - w(ip)*N(ip,:)'*-drc_dFEOH;
	
						dqC_dC(:,:,6,1) = dqC_dC(:,:,6,1) - w(ip)*N(ip,:)'*drc_dH;
						dqC_dC(:,:,6,5) = dqC_dC(:,:,6,5) - w(ip)*N(ip,:)'*drc_dFE;
						dqC_dC(:,:,6,6) = dqC_dC(:,:,6,6) - w(ip)*N(ip,:)'*(drc_dFEOH-drc2_dFEOH);
					end

					% lumped integration factor
					C_Lumped = C_Lumped + w(ip)*N(ip,:)';
				end


				for i=1:length(dofsE) %% lumped integrations
% 					for s=1:obj.n_species
% 						q_E(i) = q_E(i) + C_Lumped(i)*obj.z(s)*C(i,s);
% 						dqE_dC(i,i,s) = dqE_dC(i,i,s) + C_Lumped(i)*obj.z(s);
% 					end

					%capacity
					for s=1:obj.n_species
						q_C(i,s) = q_C(i,s) +  C_Lumped(i)/dt*(C(i,s)-C_OLD(i,s));
						dqC_dC(i,i,s,s) = dqC_dC(i,i,s,s) + C_Lumped(i)/dt;
					end

					
					if (obj.Lumped(1)) %water
						rwc = obj.k(1);
						kdum = 1e10;
						phlim=1e-16;

						rw = rwc*(1e-8-(max(phlim,C(i,1)))*(max(phlim,C(i,2))));
						drw_dH  = rwc*(-(max(phlim,C(i,2))));
						drw_dOH = rwc*(-(max(phlim,C(i,1))));

						q_C(i,1) = q_C(i,1) - C_Lumped(i)*rw;
						q_C(i,2) = q_C(i,2) - C_Lumped(i)*rw;
	
						dqC_dC(i,i,1,1) = dqC_dC(i,i,1,1) - C_Lumped(i)*drw_dH;
						dqC_dC(i,i,1,2) = dqC_dC(i,i,1,2) - C_Lumped(i)*drw_dOH;
						dqC_dC(i,i,2,1) = dqC_dC(i,i,2,1) - C_Lumped(i)*drw_dH;
						dqC_dC(i,i,2,2) = dqC_dC(i,i,2,2) - C_Lumped(i)*drw_dOH;
					end
					if (obj.Lumped(2)) %Fe
						rcc = [obj.k(2), obj.k(3), obj.k(4)];
						rc = rcc(1)*C(i,5)-rcc(2)*C(i,6)*C(i,1);
						drc_dH = -rcc(2)*C(i,6);
						drc_dFE = rcc(1);
						drc_dFEOH = -rcc(2)*C(i,1);
	
						rc2 = rcc(3)*C(i,6);
						drc2_dFEOH = rcc(3);

						q_C(i,1) = q_C(i,1) - C_Lumped(i)*(rc+rc2);
						q_C(i,5) = q_C(i,5) - C_Lumped(i)*-rc;
						q_C(i,6) = q_C(i,6) - C_Lumped(i)*(rc-rc2);
	
						dqC_dC(i,i,1,1) = dqC_dC(i,i,1,1) - C_Lumped(i)*drc_dH;
						dqC_dC(i,i,1,5) = dqC_dC(i,i,1,5) - C_Lumped(i)*drc_dFE;
						dqC_dC(i,i,1,6) = dqC_dC(i,i,1,6) - C_Lumped(i)*(drc_dFEOH+drc2_dFEOH);
	
						dqC_dC(i,i,5,1) = dqC_dC(i,i,5,1) - C_Lumped(i)*-drc_dH;
						dqC_dC(i,i,5,5) = dqC_dC(i,i,5,5) - C_Lumped(i)*-drc_dFE;
						dqC_dC(i,i,5,6) = dqC_dC(i,i,5,6) - C_Lumped(i)*-drc_dFEOH;
	
						dqC_dC(i,i,6,1) = dqC_dC(i,i,6,1) - C_Lumped(i)*drc_dH;
						dqC_dC(i,i,6,5) = dqC_dC(i,i,6,5) - C_Lumped(i)*drc_dFE;
						dqC_dC(i,i,6,6) = dqC_dC(i,i,6,6) - C_Lumped(i)*(drc_dFEOH-drc2_dFEOH);
					end
				end

				%% add to sparse allocation vectors
				for s1=1:obj.n_species
					for s2=1:obj.n_species
						[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),dofsC(:,s2));
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = dqC_dC(:,:,s1,s2);
						kmat = [kmat; tmp(:)];
					end

					[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),dofsE);
					dofmatX = [dofmatX; dofmatxloc(:)];
					dofmatY = [dofmatY; dofmatyloc(:)];
					tmp = squeeze(dqC_dE(:,:,s1));
					kmat = [kmat; tmp(:)];

					[dofmatxloc,dofmatyloc] = ndgrid(dofsE,dofsC(:,s1));
					dofmatX = [dofmatX; dofmatxloc(:)];
					dofmatY = [dofmatY; dofmatyloc(:)];
					tmp = squeeze(dqE_dC(:,:,s1));
					kmat = [kmat; tmp(:)];

					fvec = [fvec; q_C(:,s1)];
					dofvec = [dofvec; dofsC(:,s1)];
				end

				[dofmatxloc,dofmatyloc] = ndgrid(dofsE,dofsE);
                dofmatX = [dofmatX; dofmatxloc(:)];
                dofmatY = [dofmatY; dofmatyloc(:)];
                kmat = [kmat; dqE_dE(:)];

				fvec = [fvec; q_E];
                dofvec = [dofvec; dofsE];
            end 
			
			%% add to stiffness matrix and force vector
            physics.fint = physics.fint + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint), 1);
            physics.K = physics.K + sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));

            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
		end

		function plotFields(obj, physics) 
			%function which plots concentrations and pH

			for el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(el,:);
				Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1+s), elnodes);
				end

                order = [1 2 3]; % [1 4 2 5 3 6]
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec(Cdofs(order,1));
				OH(el,:) = physics.StateVec(Cdofs(order,2));
				FE(el,:) = physics.StateVec(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec(Cdofs(order,6));
				O(el,:) = physics.StateVec(Cdofs(order,7));
				E(el,:) = physics.StateVec(Edofs(order));
			end

			subplot(2,2,3)
			patch(X',Y',O',O','FaceColor','interp','EdgeColor','Interp');
			title("$C_{O2}$",'Interpreter','latex')
			colorbar

			subplot(2,2,4)
			patch(X',Y',E',E','FaceColor','interp','EdgeColor','Interp');
			title("$\varphi$",'Interpreter','latex')
			colorbar

			subplot(2,2,1)
			H(H<=0)=NaN;
			H(~isnan(H)) = -log10(H(~isnan(H))/1000); %
			patch(X',Y',H',H','FaceColor','interp','EdgeColor','Interp');
			title("pH",'Interpreter','latex')
			colorbar
			%caxis([2 12])

			subplot(2,2,2)
			patch(X',Y',FEOH',FEOH','FaceColor','interp','EdgeColor','Interp');
			title("$C_{FeOH}$",'Interpreter','latex')
			colorbar
			
		end
        
		function plotFieldsSurf(obj, physics)  %function which plots concentrations
			for el=1:size(obj.mesh.Elementgroups{obj.myGroupIndex}.Elems, 1)
                elnodes =physics.mesh.Elementgroups{obj.myGroupIndex}.Elems(el,:);
				Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
				for s=1:obj.n_species
					Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1+s), elnodes);
				end

                order = [1 2 3]; % [1 4 2 5 3 6]
                X(el,:) = physics.mesh.Nodes(elnodes(order),1);
                Y(el,:) = physics.mesh.Nodes(elnodes(order),2);
                H(el,:) = physics.StateVec(Cdofs(order,1));
				OH(el,:) = physics.StateVec(Cdofs(order,2));%
				FE(el,:) = physics.StateVec(Cdofs(order,5));
				FEOH(el,:) = physics.StateVec(Cdofs(order,6));
				O(el,:) = physics.StateVec(Cdofs(order,7));
				E(el,:) = physics.StateVec(Edofs(order));
			end
            %patch(X',Y',Z','EdgeColor','None','FaceColor','interp');
			subplot(2,2,1)
			patch(X',Y',H',H','FaceColor','interp');
			title("C_{H^+}")

			subplot(2,2,2)
			patch(X',Y',OH',OH','FaceColor','interp');
			title("C_{OH^-}")
				
			subplot(2,2,3)
			patch(X',Y',FE',FE','FaceColor','interp');
			title("FeOH^+")

			subplot(2,2,4)
			patch(X',Y',O',O','FaceColor','interp');
			title("O_2")
		end
        
    end
end

