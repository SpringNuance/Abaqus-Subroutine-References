classdef ElectrolyteInterface < BaseModel
    %Implements hydrogen and oxygen reactions at a cathodic,
	% and corrosion reactions at an anodic surface, also has the ability to
	% enforce charge conservation. Input example:
	%physics_in{2}.type = "ElectrolyteInterface";
	%physics_in{2}.Anode = "Anode";
	%physics_in{2}.Cathode = "Cathode";
	%physics_in{2}.k = [	1e-1/F_const,	1e-1/F_const,	0.5,	-0.4;  % Fe <-> Fe2+ (anode)
	%					1e-4/F_const,	1e-6/F_const,	0.5,	0;     % H+ <->H2    (cathode)
	%					1e-6/F_const,	1e-6/F_const,	0.5,	0.4;   % O2 <->OH-   (cathode)
	%					]; %reaction constants k, k', alpha, E_eq
	%physics_in{2}.ChargeConserve = true; %Enforce charge-conservation conditions
	%physics_in{2}.Em = 0;	%Metal potential (unused if chargeconserve=true
	%physics_in{2}.Lumped = [1 1 1];	%Use lumped integration to stabilise surface reactions?
    
    properties
        mesh			%pointer to mesh
        myName			%name of this model
        AnodeGroup		%name of the element group associated with this model
        AnodeGroupIndex	%index of element group
		CathodeGroup		%name of the element group associated with this model
        CathodeGroupIndex	%index of element group
        dofSpace		%pointer to degrees of freedom
        dofTypeIndices %vector containing the dof numbers associated with this model
		EmTypeIndex    %index for the metal potential

		Em		%metal potential if non-charge conserving
		ChargeConserve		%Charge conserving or imposed potential
		Lumped	%flags indicating whether to use lumped integration

		F_const = 96485.3329;		%faraday constant
		R_const = 8.31446261815324;	%gas constant
		T_const = 293.15;			%temperature

		I_anode;		%total anodic current
		I_Cathode1;		%total cathodic current
		I_Cathode2;		%total cathodic current

		k
		n_species					%number of ion species
    end
    
    methods
        function obj = ElectrolyteInterface(mesh, physics, inputs)
            %% save inputs to object
            obj.myName = "ElectrolyteInterface";
            disp("Initializing "+obj.myName)
            obj.mesh = mesh;
            obj.AnodeGroup = inputs.Anode;
            obj.AnodeGroupIndex = obj.mesh.getGroupIndex(obj.AnodeGroup);
            obj.CathodeGroup = inputs.Cathode;
            obj.CathodeGroupIndex = obj.mesh.getGroupIndex(obj.CathodeGroup);

            obj.dofSpace = physics.dofSpace;
            
            %% create relevant dofs	                       %1     2     3    4    5     6      7    8
            obj.dofTypeIndices = obj.dofSpace.addDofType({"Epot","H", "OH","Na", "Cl","Fe","FeOH","O2"});
            obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.AnodeGroupIndex));
			obj.dofSpace.addDofs(obj.dofTypeIndices, obj.mesh.GetAllNodesForGroup(obj.CathodeGroupIndex));
			obj.EmTypeIndex = obj.dofSpace.addDofType({"Em"});
			obj.dofSpace.addDofs(obj.EmTypeIndex, 1);
            
            %% get parameters
			obj.Lumped = inputs.Lumped;
			obj.k = inputs.k;

			obj.ChargeConserve = inputs.ChargeConserve;
			obj.Em = inputs.Em;

			obj.n_species = 7;
			obj.I_anode = 0;
			obj.I_Cathode1 = 0;
			obj.I_Cathode2 = 0;
        end
        
        function getKf(obj, physics)
			% Force vector and tangential matrix assembly procedure

            fprintf("        ElectrolyteInterface get Matrix:")
            t = tic;
            
            dt = physics.dt;

			EmDof = obj.dofSpace.getDofIndices(obj.EmTypeIndex, 1);
			I_an = 0;
			I_cat1 = 0;
			I_cat2 = 0;
			if (obj.ChargeConserve)
				EM = physics.StateVec(EmDof);
				obj.Em = EM;
			else
				physics.condofs = [physics.condofs; EmDof];
                physics.convals = [physics.convals; obj.Em];
				EM = obj.Em;
			end

			%stiffness matrix
            dofmatX = [];
            dofmatY = [];
            kmat = [];
            fvec = [];
            dofvec = [];

			Svec = physics.StateVec;
            SvecOld = physics.StateVec_Old;

			for AC=1:2 %perform for both surfaces
				if AC==1
					Group = obj.AnodeGroupIndex;
					surface = "Anode";
				else
					Group = obj.CathodeGroupIndex;
					surface = "Cathode";
				end

				parfor n_el=1:size(obj.mesh.Elementgroups{Group}.Elems, 1)
					%get nodes and element shape functions
                	Elem_Nodes = obj.mesh.getNodes(Group, n_el); 
                	[N, ~, w] = obj.mesh.getVals(Group, n_el);
					ipcoords = obj.mesh.getIPCoords(Group, n_el);
					w = w.*2.*pi.*ipcoords(1,:);
	
					%get relevant dofs and nodal values
					dofsE = obj.dofSpace.getDofIndices(obj.dofTypeIndices(1), Elem_Nodes);
					dofsC = zeros(length(dofsE), obj.n_species);
					for s=1:obj.n_species
						dofsC(:,s) = obj.dofSpace.getDofIndices(obj.dofTypeIndices(s+1), Elem_Nodes);
					end
	
					C = Svec(dofsC);
					E = Svec(dofsE);
	
					C_OLD = SvecOld(dofsC);
					E_OLD = SvecOld(dofsE);
	
					%initialize arrays
					q_C    = zeros(length(dofsE), obj.n_species);
					dqC_dC = zeros(length(dofsE), length(dofsE), obj.n_species, obj.n_species);
					dqC_dE = zeros(length(dofsE), length(dofsE), obj.n_species);
					dqC_dEM = zeros(length(dofsE), 1, obj.n_species);
	
					q_Em = 0;
					dqEm_dEm = zeros(1, 1);
					dqEm_dE = zeros(1, length(dofsE));
					dqEm_dC = zeros(1, length(dofsE), obj.n_species);
	
					C_Lumped = zeros(length(dofsE), 1);
					for ip=1:length(w)
	
						%% surface reactions
						CH = N(ip,:)*C(:,1);
						COH = N(ip,:)*C(:,2);
						CFE = N(ip,:)*C(:,5);
						CFEOH = N(ip,:)*C(:,6);
						CO2 = N(ip,:)*C(:,7);
						phil = N(ip,:)*E;
	
						[react, dreact, products, ions] = obj.reactions(CH, COH, CFE, CFEOH, CO2, EM, phil, surface);
	
						% total reactions
						for r=1:length(ions)
							for s=1:obj.n_species
								q_C(:,s) = q_C(:,s) - w(ip)*N(ip,:)'*(react(r,1)-react(r,2))*products(r,s) *(1-obj.Lumped(r));
								for n=1:obj.n_species
									dqC_dC(:,:,s,n) = dqC_dC(:,:,s,n)-w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,2+n)-dreact(r,2,2+n))*products(r,s)*(1-obj.Lumped(r));
								end
								dqC_dE(:,:,s)   = dqC_dE(:,:,s)  - w(ip)*N(ip,:)'*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*products(r,s)*(1-obj.Lumped(r));
								dqC_dEM(:,1,s)  = dqC_dEM(:,1,s) - w(ip)*N(ip,:)'*(dreact(r,1,2)-dreact(r,2,2))*products(r,s)*(1-obj.Lumped(r));
							end
							q_Em      = q_Em    + w(ip)*(react(r,1)-react(r,2))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
							dqEm_dEm = dqEm_dEm + w(ip)*(dreact(r,1,2)-dreact(r,2,2))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
							dqEm_dE = dqEm_dE   + w(ip)*N(ip,:)*(dreact(r,1,1)-dreact(r,2,1))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
							for n=1:obj.n_species
								dqEm_dC(1,:,n) = dqEm_dC(1,:,n)+ w(ip)*N(ip,:)*(dreact(r,1,2+n)-dreact(r,2,2+n))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
							end
							if (surface == "Anode")
								I_an = I_an + w(ip)*(react(r,1)-react(r,2))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
							else
								if r==2
								I_cat1 = I_cat1 + w(ip)*(react(r,1)-react(r,2))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
								else
								I_cat2 = I_cat2 + w(ip)*(react(r,1)-react(r,2))*ions(r)*(1-obj.Lumped(r))*obj.F_const;
								end
							end
						end
	
						% lumped integration weight
						C_Lumped = C_Lumped + w(ip)*N(ip,:)';
					end
	
	
					for i=1:length(dofsE) %% lumped integrations
						[react, dreact, products, ions] = obj.reactions(C(i,1), C(i,2),C(i,5),C(i,6),C(i,7), EM, E(i), surface);
	
						for r=1:length(ions)
							for s=1:obj.n_species
								q_C(i,s) = q_C(i,s) - C_Lumped(i)*(react(r,1)-react(r,2))*products(r,s) *obj.Lumped(r);
								for n=1:obj.n_species
									dqC_dC(i,i,s,n) = dqC_dC(i,i,s,n)-C_Lumped(i)*(dreact(r,1,2+n)-dreact(r,2,2+n))*products(r,s)*obj.Lumped(r);
								end
								dqC_dE(i,i,s)   = dqC_dE(i,i,s)  - C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*products(r,s)*obj.Lumped(r);
								dqC_dEM(i,1,s)  = dqC_dEM(i,1,s) - C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*products(r,s)*obj.Lumped(r);
							end
							q_Em      = q_Em    + C_Lumped(i)*(react(r,1)-react(r,2))*ions(r)*obj.Lumped(r)*obj.F_const;
							dqEm_dEm = dqEm_dEm + C_Lumped(i)*(dreact(r,1,2)-dreact(r,2,2))*ions(r)*obj.Lumped(r)*obj.F_const;
							dqEm_dE(i) = dqEm_dE(i)   + C_Lumped(i)*(dreact(r,1,1)-dreact(r,2,1))*ions(r)*obj.Lumped(r)*obj.F_const;
							for n=1:obj.n_species
								dqEm_dC(1,i,n) = dqEm_dC(1,i,n)+ C_Lumped(i)*(dreact(r,1,2+n)-dreact(r,2,2+n))*ions(r)*obj.Lumped(r)*obj.F_const;
							end
							if (surface == "Anode")
								I_an = I_an + C_Lumped(i)*(react(r,1)-react(r,2))*ions(r)*obj.Lumped(r)*obj.F_const;
							else
								if r==2
									I_cat1 = I_cat1 + C_Lumped(i)*(react(r,1)-react(r,2))*ions(r)*obj.Lumped(r)*obj.F_const;
								else
									I_cat2 = I_cat2 + C_Lumped(i)*(react(r,1)-react(r,2))*ions(r)*obj.Lumped(r)*obj.F_const;
								end
							end
						end
					end
	
					%% save to sparse allocation vectors
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
						tmp = dqC_dE(:,:,s1);
						kmat = [kmat; tmp(:)];
	
						[dofmatxloc,dofmatyloc] = ndgrid(dofsC(:,s1),EmDof);
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = dqC_dEM(:,:,s1);
						kmat = [kmat; tmp(:)];
	
						[dofmatxloc,dofmatyloc] = ndgrid(EmDof,dofsC(:,s1));
						dofmatX = [dofmatX; dofmatxloc(:)];
						dofmatY = [dofmatY; dofmatyloc(:)];
						tmp = dqEm_dC(:,:,s1);
						kmat = [kmat; tmp(:)];
	
						fvec = [fvec; q_C(:,s1)];
						dofvec = [dofvec; dofsC(:,s1)];
					end
	
					[dofmatxloc,dofmatyloc] = ndgrid(EmDof,EmDof);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqEm_dEm(:)];
	
					[dofmatxloc,dofmatyloc] = ndgrid(EmDof,dofsE);
                	dofmatX = [dofmatX; dofmatxloc(:)];
                	dofmatY = [dofmatY; dofmatyloc(:)];
                	kmat = [kmat; dqEm_dE(:)];
	
					fvec = [fvec; q_Em];
                	dofvec = [dofvec; EmDof];
				end 
			end
			obj.I_anode   = I_an;
			obj.I_Cathode1 = I_cat1;
			obj.I_Cathode2 = I_cat2;

			% add to stiffness matrix and force vector
            physics.fint = physics.fint + sparse(dofvec, 0*dofvec+1, fvec, length(physics.fint), 1);
            physics.K = physics.K + sparse(dofmatX, dofmatY, kmat, length(physics.fint),length(physics.fint));

            tElapsed = toc(t);
            fprintf("            (Assemble time:"+string(tElapsed)+")\n");
			fprintf("            Em: " +string(obj.Em)+ "  \n            Anodic Current: " +string(obj.I_anode)+ "  \n            Cathodic Current: "+string(obj.I_Cathode1+obj.I_Cathode2) +"\n            Mismatch: "+string(obj.I_anode+obj.I_Cathode1++obj.I_Cathode2)+"\n")
		end

		function [rout1,rout2,rout3,xout,Eout] = ReactionsVsX(obj, physics)
			%returns the local reaction rates at the surface
			xout = [];
			rout1 = [];
			rout2 = [];
			rout3 = [];
			Eout = [];
			cnt=0;
			for AC=1:2
				if AC==1
					Group = obj.AnodeGroupIndex;
					surface = "Anode";
				else
					Group = obj.CathodeGroupIndex;
					surface = "Cathode";
				end
		
				for el=1:size(obj.mesh.Elementgroups{Group}.Elems, 1)
					cnt=cnt+1;
                	elnodes =physics.mesh.Elementgroups{Group}.Elems(el,:);
					Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
					for s=1:obj.n_species
						Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1+s), elnodes);
					end
	
                	order = [1 2 3];
                	X(cnt,:) = [physics.mesh.Nodes(elnodes(order),1)];
                	Y(cnt,:) = [physics.mesh.Nodes(elnodes(order),2)];	

                	H(cnt,:) = physics.StateVec(Cdofs(order,1));
					OH(cnt,:) = physics.StateVec(Cdofs(order,2));%
					FE(cnt,:) = physics.StateVec(Cdofs(order,5));
					FEOH(cnt,:) = physics.StateVec(Cdofs(order,6));
					O(cnt,:) = physics.StateVec(Cdofs(order,7));
					E(cnt,:) = physics.StateVec(Edofs(order));
	
					EmDof = obj.dofSpace.getDofIndices(obj.EmTypeIndex, 1);
					EM = physics.StateVec(EmDof);
	
					for i=1:length(order)
						[react, ~, ~, ions] = obj.reactions(H(cnt,i), OH(cnt,i),FE(cnt,i),FEOH(cnt,i),O(cnt,i) ...
						                     	             , EM, E(cnt,i), surface);
						for j=1:3
							r(cnt,i,j) = (react(j,1)-react(j,2))*ions(j);
						end
					end
					Eplot(cnt,:) = EM-E(cnt,:);

					xout = [xout; X(cnt,:)];
					Eout = [Eout; Eplot(cnt,:)];
					rout1= [rout1; r(cnt,:,1)];
					rout2= [rout2; r(cnt,:,2)];
					rout3= [rout3; r(cnt,:,3)];
				end
			end

			[~, idx, ~] = unique(xout);
			xout = xout(idx);
			Eout = Eout(idx);
			rout1 = rout1(idx);
			rout2 = rout2(idx);
			rout3 = rout3(idx);

			[~,idx] = sort(xout);
			xout = xout(idx);
			Eout = Eout(idx);
			rout1 = rout1(idx)*2*obj.F_const;
			rout2 = rout2(idx)*2*obj.F_const;
			rout3 = rout3(idx)*4*obj.F_const;
		end

		function plotReactions(obj, physics) %plots surface reactions
			%plots reaction rates for individual reactions (plotting
			%performed based on integration-point reaction rates)

			cnt=0;
			for AC=1:2
				if AC==1
					Group = obj.AnodeGroupIndex;
					surface = "Anode";
				else
					Group = obj.CathodeGroupIndex;
					surface = "Cathode";
				end
		
				for el=1:size(obj.mesh.Elementgroups{Group}.Elems, 1)
					cnt=cnt+1;
                	elnodes =physics.mesh.Elementgroups{Group}.Elems(el,:);
					Edofs = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1), elnodes);
					for s=1:obj.n_species
						Cdofs(:,s) = physics.dofSpace.getDofIndices(obj.dofTypeIndices(1+s), elnodes);
					end
	
                	order = [1 2 3];
                	X(cnt,:) = [physics.mesh.Nodes(elnodes(order),1);NaN];
                	Y(cnt,:) = [physics.mesh.Nodes(elnodes(order),2);NaN];	

                	H(cnt,:) = physics.StateVec(Cdofs(order,1));
					OH(cnt,:) = physics.StateVec(Cdofs(order,2));%
					FE(cnt,:) = physics.StateVec(Cdofs(order,5));
					FEOH(cnt,:) = physics.StateVec(Cdofs(order,6));
					O(cnt,:) = physics.StateVec(Cdofs(order,7));
					E(cnt,:) = physics.StateVec(Edofs(order));
	
					EmDof = obj.dofSpace.getDofIndices(obj.EmTypeIndex, 1);
					EM = physics.StateVec(EmDof);
	
					for i=1:length(order)
						[react, ~, ~, ions] = obj.reactions(H(cnt,i), OH(cnt,i),FE(cnt,i),FEOH(cnt,i),O(cnt,i) ...
						                     	             , EM, E(cnt,i), surface);
						for j=1:3
							r(cnt,i,j) = (react(j,1)-react(j,2))*ions(j);
						end
					end
					r(cnt,4,:) = NaN;
					Eplot(cnt,1:3) = EM-E(cnt,:);
					Eplot(cnt,4) = NaN;
				end
			end
            %Iron
			subplot(2,2,1)
			fill3(X',Y',r(:,:,1)',r(:,:,1)','FaceColor','interp','EdgeColor','interp','LineWidth',3);
			title("\nu_{Fe}")
			view(2)
			hold on
			colorbar

			%Hydrogen
			subplot(2,2,2)
			fill3(X',Y',r(:,:,2)',r(:,:,2)','FaceColor','interp','EdgeColor','interp','LineWidth',3);
			title("\nu_{H}")
			view(2)
			hold on
			colorbar
				
			%Oxygen
			subplot(2,2,3)
			fill3(X',Y',r(:,:,3)',r(:,:,3)','FaceColor','interp','EdgeColor','interp','LineWidth',3);
			title("\nu_{O2}")
			view(2)
			hold on
			colorbar

			%Overpotential
			subplot(2,2,4)
			fill3(X',Y',Eplot',Eplot','FaceColor','interp','EdgeColor','interp','LineWidth',3);
			title("$E_m-\varphi$",'Interpreter','latex')
			view(2)
			hold on
			colorbar
		end
        
		function [react, dreact, products, ions] = reactions(obj, CH, COH, CFE, CFEOH, CO2, EM, phil, surface)
			%function calculating the reactions based on input values (either nodal
			%or integration point), outputs the reaction rate, derivative of
			%reaction rate, and count of ions involved with reaction
			if (surface == "Anode")
					products(1,:) = [0, 0, 0, 0, -1, 0, 0];     %Fe2+ +2e-   <-> Fe 
					products(2,:) = [0, 0, 0, 0, 0, 0, 0];      %2H+ +2e-    <-> H2 
					products(3,:) = [0, 0, 0, 0, 0, 0, 0];      %O2+2H2O+4e- <-> 4OH-

					ions = [2 0 0];
			else
					products(1,:) = [0, 0, 0, 0, 0, 0, 0];      %Fe2+ +2e-   <-> Fe 
					products(2,:) = [-2, 0, 0, 0, 0, 0, 0];      %2H+ +2e-    <-> H2 
					products(3,:) = [0, 4, 0, 0, 0, 0, -1];      %O2+2H2O+4e- <-> 4OH-

					ions = [0 2 4];
			end

			react = zeros(3,2); %reaction, forward/backward
			dreact= zeros(3,2,2+obj.n_species); %reaction, forward/backwards, E/Em/species

			%corrosion
			react(1,1)    = obj.k(1,1)*CFE*exp(-obj.k(1,3)*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(1,1,1) = obj.k(1,1)*CFE*exp(-obj.k(1,3)*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(1,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(1,1,2) = obj.k(1,1)*CFE*exp(-obj.k(1,3)*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(1,3)*(1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(1,1,7) = obj.k(1,1)*exp(-obj.k(1,3)*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
			
			react(1,2)    = obj.k(1,2)*exp((1-obj.k(1,3))*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(1,2,1) = obj.k(1,2)*exp((1-obj.k(1,3))*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(1,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(1,2,2) = obj.k(1,2)*exp((1-obj.k(1,3))*(EM-phil-obj.k(1,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(1,3))*(1)*obj.F_const/obj.R_const/obj.T_const);

			%H2 (no backwards assuming H2 disappears)
			react(2,1)    = obj.k(2,1)*CH*exp(-obj.k(2,3)*(EM-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(2,1,3) = obj.k(2,1)*exp(-obj.k(2,3)*(EM-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(2,1,1) = obj.k(2,1)*CH *exp(-obj.k(2,3)*(EM-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(2,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(2,1,2) = obj.k(2,1)*CH *exp(-obj.k(2,3)*(EM-phil-obj.k(2,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(2,3)*(1)*obj.F_const/obj.R_const/obj.T_const);

			%O2
			react(3,1)    = obj.k(3,1)*CO2*exp(-obj.k(3,3)*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,1,9) = obj.k(3,1)*exp(-obj.k(3,3)*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,1,1) = obj.k(3,1)*CO2 *exp(-obj.k(3,3)*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(3,3)*(-1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,1,2) = obj.k(3,1)*CO2 *exp(-obj.k(3,3)*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const)*(-obj.k(3,3)*(1)*obj.F_const/obj.R_const/obj.T_const);

			react(3,2)    = obj.k(3,2)*COH*exp((1-obj.k(3,3))*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,2,4) = obj.k(3,2)*exp((1-obj.k(3,3))*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,2,1) = obj.k(3,2)*COH *exp((1-obj.k(3,3))*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(3,3))*(-1)*obj.F_const/obj.R_const/obj.T_const);
			dreact(3,2,2) = obj.k(3,2)*COH *exp((1-obj.k(3,3))*(EM-phil-obj.k(3,4))*obj.F_const/obj.R_const/obj.T_const)*((1-obj.k(3,3))*(1)*obj.F_const/obj.R_const/obj.T_const);
		end
        
    end
end

