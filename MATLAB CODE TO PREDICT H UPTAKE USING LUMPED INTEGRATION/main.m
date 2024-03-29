% Finite element simulation code accompanying T. Hageman, E. Martinez-Paneda, 
% Stabilising Effects of Lumped Integration Schemes for the Simulation of 
% Metal-Electrolyte Reactions. Journal of The Electrochemical Society (2023), 
% http://www.doi.org/10.1149/1945-7111/acb971
%
% This performs the simulation of a metal domain interacting with an electrolyte, 
% and demonstrates the beneficial effects of lumped integration. If this
% code is used, either directly or after modification, please cite 
% http://www.doi.org/10.1149/1945-7111/acb971 in any resulting research. 
%
% Usage: Set the parameters for the case being considered within this file,
% main.m, and run this file. Other functions contained within this code are
% automatically called, and outputs are automatically saved to the defined
% save folder. Verified to work with matlab 2021b and 2022a
%
% DISCLAIMER: While this code has been cross-verified with comparison to
% simulation results from the commercial finite element package COMSOL, it
% can not be guaranteed to be error-free. Before using this code for
% relevant or critical applications, especially when simulating cases not 
% directly included, please perform your own verification. The authors are 
% not repsonsible for any issues arrising from mistakes within this matlab code.
close all; clear all; clc
addpath(genpath('./Models'))
addpath(genpath('./Shapes'))

maxNumCompThreads(8); %Number of cpu cores to use during simulation (code uses parfor parrallellisation, thus can only run on a single node)
delete(gcp('nocreate'))
parpool('threads')
%% Input parameters used in sweeps
	Em = 0;	%metal potential
	fprintf('Starting job:'+string(Em)+'\n');

	k = [1e-4,	1e-10,	0.5,	0;   %reaction constants, [k, k', alpha, E_eq] for	Acidic Volmer
	   	1e-10,	0,		0.3,	0;	 %											Acidic Heyrovsky
	   	1e-6,	0,		0,		0;	 %											Tafel
	   	1e3,	7e7,	0,		0;	 %											Absorbtion
	   	1e-8,	1e-13,	0.5,	0;	 %											Basic Volmer
	   	1e-10,	1e-14,	0.3,	0;   %											Basic Heyrovsky
	   	3e-5/(2*96485.3329),3e-5/(2*96485.3329), 0.5, -0.4];   %				Corrosion
	sname = "Em_"+string(Em);

	savefolder = "./Results/"+sname; % Folder in which to save output files
	mkdir(savefolder);
	savefolder=savefolder+"/";

%% input properties used, sorted per physical model
	restart = false;
	restart_num = -1;
	if restart == false
    	%% mesh properties
		% Finer control over the geometry and element size is possible
		% within ./@Mesh/Geometry_Generator.m 
    	mesh_in.Lx    = 10e-3;	%Domain length [m]
    	mesh_in.Ly    = 10e-3;	%Domain Height [m]
    	mesh_in.ipcount1D = 3;	%numer of integration points per dimension
    	mesh_in.zeroWeight = true; %Include zero-weight integration points at element boundaries (for post-processing)
	
    	%% physics models
		% Linear elastic material description for metal domain
    	physics_in{1}.type = "LinearElastic";	
    	physics_in{1}.Egroup = "Metal";	
    	physics_in{1}.young = 200e9;	% Youngs modulus [Pa]
    	physics_in{1}.poisson = 0.3;	% Poisson ratio [-]
	
		% Hydrogen diffusion within the metal domain
    	physics_in{2}.type = "HydrogenDiffusion";
		physics_in{2}.Egroup = "Metal";
    	physics_in{2}.DL = 1e-9;	% Lattice difusivity
		physics_in{2}.NL = 1e6;		% Amount of interstitial lattice sites [mol/m^3]
	
		% Displacement constraints for metal domain
    	physics_in{3}.type = "Constrainer";
    	physics_in{3}.Egroup = "M_Bottom";
    	physics_in{3}.dofs = {"dx"};
    	physics_in{3}.conVal = [0];
	
    	physics_in{4}.type = "Constrainer";
    	physics_in{4}.Egroup = "M_Bottom";
    	physics_in{4}.dofs = {"dy"};
    	physics_in{4}.conVal = [0];
	
		physics_in{5}.type = "Constrainer";
    	physics_in{5}.Egroup = "M_Top";
    	physics_in{5}.dofs = {"dx"};
    	physics_in{5}.conVal = [0];
	
    	physics_in{6}.type = "Constrainer";
    	physics_in{6}.Egroup = "M_Top";
    	physics_in{6}.dofs = {"dy"};
    	physics_in{6}.conVal = [0.01e-3];
	
		% Nernst-planck, electroneutrality, and volume reactions for electrolyte
		physics_in{7}.type = "Electrolyte";
    	physics_in{7}.Egroup = "Electrolyte";
    	physics_in{7}.D = [9.3; 5.3; 1.3; 2; 1.4; 1]*1e-9;  % Diffusion coefficients [m/s] for ions: H OH Na Cl Fe FeOH
		physics_in{7}.z = [1; -1; 1; -1; 2; 1];	% ionic charges
		physics_in{7}.pH0 = 5;	% Initial condition pH
		physics_in{7}.NaCl = 0.6e3; % Initial concentration of NaCl
		physics_in{7}.Lumped = [true; true]; % Flag for using lumped integration for water auto-ionisation and metal-ion reactions
		physics_in{7}.k = [1e6; 1e-1; 1e-3; 1e-3]; % Reaction constants k_eq, k_fe, k_fe', k_feoh
	
		% Boundary constraints for electrolyte
		initH = 1000*10^(-physics_in{7}.pH0);
		initOH = 1000*10^(-14+physics_in{7}.pH0);
		initCl = physics_in{7}.NaCl;
		initNa = initCl-initH+initOH;
	
		physics_in{8}.type = "Constrainer";
    	physics_in{8}.Egroup = "E_Left";
    	physics_in{8}.dofs = {"Epot";"H";"OH";"Na";"Cl";"Fe";"FeOH"};
    	physics_in{8}.conVal = [0; initH; initOH; initNa; initCl; 0; 0];
	
		% Metal-electrolyte interface
		physics_in{9}.type = "ElectrolyteInterface";
    	physics_in{9}.Egroup = "Interface";
		physics_in{9}.NAds = 1e-3; % Concentration of surface sites [mol/m^2]
		physics_in{9}.k = k; %Reaction constants
		physics_in{9}.NL = physics_in{2}.NL; % Concentration of interstitial lattice sites [mol/m^3]
		physics_in{9}.Em = Em; % Metal Potential [V_SHE]
		physics_in{9}.Lumped = [1 1 1 1 1 1 1]; %Flags to enable lumped integration on a per-reaction basis
	
	
    	%% solver inputs
    	solver_in.maxIt = 100;	% Maximum number of Newton-Raphson iterations per time step
    	solver_in.Conv = 1e-9;	% Convergence criterion (relative energy based)
    	solver_in.tiny = 1e-9;	% Convergence criterion (absolute energy based)
    	solver_in.linesearch = true; % Flag to allow for a linear line-search
    	solver_in.linesearchLims = [0.1 1]; % Bounding box for line-search algorith
	
    	%% initialization

		% Generate mesh
    	mesh = Mesh(mesh_in);
    	%mesh.plot(true, true, true);  %(enable to plot the used mesh at the start)
    	mesh.check();
	
		%Initialize physics models
    	physics = Physics(mesh, physics_in);
	
		% parameters for time loop
    	dt = 30; %initial time increment
    	physics.time = 0;
	
		tmax = 60*60*24*365*50; % total time to simulation
    	n_max = 1000*24*360; % maximum number of time steps if tmax is not reached
    	solver = Solver(physics, solver_in);

		%Vectors to save time-dependent behaviour into
    	tvec = 0;
		CL_vec = 0;
		Cmax_vec = 0;
	
    	startstep = 1;
	else
		% Operators to restart from a recovery file
    	filename = savefolder+string(restart_num);
    	load(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");
    	startstep = restart_num+1;
		solver.linesearch = false;
		solver.linesearchLims = [0.1 1];
	end
	
%% perform time-dependent simulatio
	for tstep = startstep:n_max
		disp("Step: "+string(tstep));
		disp("Time: "+string(physics.time));
		physics.dt = dt*1.05^(tstep-1);
		disp("dTime: "+string(physics.dt));
        
		% solve current time increment
		solver.Solve(); 
        	
		% update results
        physics.time = physics.time+physics.dt;
        tvec(end+1) = tvec(end)+physics.dt;
		CL_vec(end+1) = physics.models{2}.CL_int./mesh.Area(1);
		Cmax_vec(end+1) = physics.models{2}.CL_max;
    	
		% plot outputs
        if mod(tstep, 1) == 0
            plotres(physics, tvec, CL_vec);
		end

		%save outputs to file
        if mod(tstep, 10) == 0
            filename = savefolder+string(tstep);
            save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");
		end
	

		if (physics.time>tmax)
			break
		end
	end

	% Save outputs one last time
	filename = savefolder+"end";
    save(filename, "mesh","physics","solver","dt","tvec","CL_vec","Cmax_vec","n_max","tmax");

function plotres(physics, tvec, CL_vec)
% PLOTRES Plots results based on the current simulation state contained
% within "physics"

    figure(42)
	clf 

	%Hydrogen concentration
    subplot(2,3,1)
        physics.PlotNodal("H",-1, "Electrolyte"); 
        title("H^+")
		colorbar

	%Hydroxide concentration
	subplot(2,3,2)
        physics.PlotNodal("OH",-1, "Electrolyte");
        title("OH^-")
		colorbar

	%Chlorine concentration
	subplot(2,3,3)
        physics.PlotNodal("CL",-1, "Metal");
        title("C_L")
		colorbar

	%Interstitial lattice concentration over time
	subplot(2,3,4)
		plot(tvec/3600, CL_vec)
		xlabel('t [hours]')
		ylabel('avarage C_L [mol/m^3]')

	%Surface coverage
	subplot(2,3,5)
		physics.PlotNodal("Theta",-1, "Interface");
		title('\theta')
		colorbar

	%Iron ion concentration
	subplot(2,3,6)
        physics.PlotNodal("Fe",-1, "Electrolyte");
        title("Fe^{2+}")
		colorbar

	%Surface reaction rates
	figure(44)
		clf
		physics.models{9}.plotReactions(physics);

	%electrolyte specific fields (pH etc.)
	figure(43)
		clf
		physics.models{7}.plotFields(physics);

     drawnow();
end

