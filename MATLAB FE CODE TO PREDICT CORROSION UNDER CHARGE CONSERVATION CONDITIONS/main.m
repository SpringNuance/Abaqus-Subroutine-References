% Finite element simulation code accompanying T. Hageman, C. Andrade & 
% E. Martínez-Pañeda, Corrosion rates under charge-conservation conditions. 
% Electrochimica Acta (2023) 
%
% This performs the simulation of a pencil-electrode test, under the 
% assumption that the metal is isolated from current sources. As such,
% the anodic corrosion reaction is restricted by the electron consumption 
% of the cathodic reactions. If this code is used, either directly or after 
% modification, please cite https://doi.org/10.1016/j.electacta.2023.142624 
% in any resulting research. 
%
% Usage: Set the parameters for the case being considered within this file,
% main.m, and run this file. Other functions contained within this code are
% automatically called, and outputs are automatically saved to the defined
% save folder. Verified to work with matlab 2021b

close all; clear all; clc

%% Varied parameters
initO2 = 0.25;		%initial and boundary oxygen concentration [mol/m^3]
initNaCL = 10^2;	%Concentration Cl- imposed at boundary [mol/m^3]
initpH = 12;		%initial and boundary pH [-]

Lfrac = 0.5e-2;			%Radius of the corrosion pit [m]
Lx = Lfrac + 1000e-3/10;	%Radius of the domain [m] (remove division/10 to replicate results from paper, included here to reduce the domain simulated for the example within the documentation)
Hfrac = 2e-2;			%depth of the corrosion pit [m]
Ly = 0.5e-2;			%Height of the domain [m]

OxLim = true;		%Flag to indicate if the oxygen boundary condition is altered within the simulations

%% set folder for saving results
savefolder = "./Results";
mkdir(savefolder);
savefolder=savefolder+"/";

%% add relevant folders to the matlab path
addpath(genpath('./Models'))	%physical models
addpath(genpath('./Shapes'))	%finite element shape functions

%% Initialize pool for multi-cpu-core computation
maxNumCompThreads(8);	
delete(gcp('nocreate'));
p = parpool('threads')

%% input properties

%check whether to start anew or resume from previous simulation files
Files=dir(savefolder);
nmax = 10*(length(Files)-2)-20;
restartfrom = nmax;
if restartfrom>0
	restart = true;
	restart_num = restartfrom;
else
	restart = false;
	restart_num = 0;
end

if restart == false
	% mesh properties
	mesh_in.type = "CorrosionPit";	%Name of the mesh generator used
	mesh_in.dxmin    = 0.1e-3;		%minimum element size 
	mesh_in.dxmax    = Ly/4;		%maximum element size
	mesh_in.Lx    = Lx;				%Domain radius [m]
	mesh_in.Ly    = Ly;				%Domain height [m]
	mesh_in.Lfrac = Lfrac;			%Corrosion pit radius
	mesh_in.Hfrac = Hfrac;			%Corrosion pit depth [m]
	mesh_in.ipcount1D = 2;			%order of numerical integration scheme
	mesh_in.zeroWeight = true;		%add zero-weight integration points to facilitate post-processing
	mesh_in.SaveName = savefolder;	%folder in which to save mesh file
	mesh_in.generate = true;		%Should a new mesh be generator, or loaded from file

	%physics models

	%Electrolyte domain
	physics_in{1}.type = "Electrolyte";
	physics_in{1}.Egroup = "Electrolyte"; 
	physics_in{1}.D = [9.3; 5.3; 1.3; 2; 1.4; 1; 1]*1e-9;  %Diffusivity of H OH Na CL Fe FeOH O2 [m/s]
	physics_in{1}.z = [1; -1; 1; -1; 2; 1; 0];				%charge of each ionic species [-]
	physics_in{1}.pH0 = initpH;								%initial pH [-]
	physics_in{1}.NaCl = initNaCL;							%Initial Cl- concentration [mol/m^3]
	physics_in{1}.O2 = initO2;								%Initial oxygen concentration [mol/m^3]
	physics_in{1}.Lumped = [true; true];		%should lumped integration be used to stabilide the water, metal volume reactions
	physics_in{1}.k = [1e6; 1e-1; 1e-3; 1e-3];	%reaction rates for the volume reactions: k_eq, k_f, k_f', k_feoh

	%Reacting surfaces
	F_const = 96485.3329; %Faraday constant
	physics_in{2}.type = "ElectrolyteInterface";
	physics_in{2}.Anode = "Anode";
	physics_in{2}.Cathode = "Cathode";
	physics_in{2}.k = [	1e-1/F_const,	1e-1/F_const,	0.5,	-0.4;  % Fe <-> Fe2+ (anode)
						1e-4/F_const,	1e-6/F_const,	0.5,	0;     % H+ <->H2    (cathode)
						1e-6/F_const,	1e-6/F_const,	0.5,	0.4;   % O2 <->OH-   (cathode)
						]; %reaction constants k, k', alpha, E_eq
	physics_in{2}.ChargeConserve = true; %Enforce charge-conservation conditions
	physics_in{2}.Em = 0;	%Metal potential (unused if chargeconserve=true
	physics_in{2}.Lumped = [1 1 1];	%Use lumped integration to stabilise surface reactions?

	%boundary conditions
	initH = 1000*10^(-physics_in{1}.pH0);
	initOH = 1000*10^(-14+physics_in{1}.pH0);
	initCl = physics_in{1}.NaCl;
	initNa = initCl-initH+initOH;

	%Top boundary
	if (OxLim) %Use a model which turns off after set amount of time
		physics_in{3}.type = "OxygenLimiter";
		physics_in{3}.Egroup = "E_Top";
		physics_in{3}.dofs = {"O2"}; %Name of degree of freedom
		physics_in{3}.conVal = [initO2];	%value of oxygen constraint
		physics_in{3}.tmax = 48*3600;	%time after which constraint is removed
	else
		physics_in{3}.type = "Constrainer";
		physics_in{3}.Egroup = "E_Top";
		physics_in{3}.dofs = {"O2"};	%Name of degree of freedom
		physics_in{3}.conVal = [initO2];	%value of oxygen constraint
	end

	%outer boundary
	physics_in{4}.type = "Constrainer";
	physics_in{4}.Egroup = "E_Right";
	physics_in{4}.dofs = {"H";"OH";"Na";"Cl";"Fe";"FeOH"};
	physics_in{4}.conVal = [initH; initOH; initNa; initCl; 0; 0];

	physics_in{5}.type = "Constrainer";
	physics_in{5}.Egroup = "E_Right";
	physics_in{5}.dofs = {"Epot"};
	physics_in{5}.conVal = [0];


	%solver inputs
	solver_in.maxIt = 20;	%max number of nonlinear iterations per time increment
	solver_in.Conv = 1e-4;	%Relative convergence criterion used
	solver_in.tiny = 1e-8;	%Absolute convergence criterion used
	solver_in.linesearch = true;	%Use a linear line-search during solution procedure
	solver_in.linesearchLims = [0.1 1]; %limits to perform linesearch within

	%% initialization
	mesh = Mesh(mesh_in);
	mesh.check();
	drawnow();

	physics = Physics(mesh, physics_in);
	physics.time = 0;

	solver = Solver(physics, solver_in);

	%time stepping
	dt = 30;
	n_max = 1000*24*360;
	tmax = 60*60*120;
	startstep = 1;

	%Initialization of vectors to save temporal data to
	tvec = 0;
	I_an_vec = 0;
	I_cat_H_vec = 0;
	I_cat_O_vec = 0;
	Em_vec = 0;
else
	%restart from previously saved file
	filename = savefolder+string(restart_num);
	load(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");
	startstep = restart_num+1;
end

%% Simulation:
%% Time stepping loop
for tstep = startstep:n_max
	%print information for commencing step
	disp("Step: "+string(tstep));
	disp("Time: "+string(physics.time));
	physics.dt = min(3600,dt*1.05^(tstep-1));
	disp("dTime: "+string(physics.dt));
	
	%solve current timestep
	solver.Solve();
	
	%save results
	physics.time = physics.time+physics.dt;
	tvec(end+1) = tvec(end)+physics.dt;
	I_an_vec(end+1) = physics.models{2}.I_anode;
	I_cat_H_vec(end+1)= physics.models{2}.I_Cathode1;
	I_cat_O_vec(end+1)= physics.models{2}.I_Cathode2;
	Em_vec(end+1)   = physics.models{2}.Em;

	if mod(tstep, 1) == 0
		plotres(physics, tvec, I_an_vec, I_cat_H_vec, I_cat_O_vec, Em_vec);
	end
	if mod(tstep, 10) == 0
		filename = savefolder+string(tstep);
		save(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");
	end

	%stop simulations once maximum time reached
	if (physics.time>tmax)
		break
	end
end

%save final results
filename = savefolder+"end";
save(filename, "mesh","physics","solver","dt","tvec","I_an_vec","I_cat_H_vec","I_cat_O_vec","Em_vec","n_max","tmax");


function plotres(physics, tvec, I_an_vec, I_cat_H_vec, I_cat_O_vec, Em_vec)
% PLOTRES Plots results based on the current simulation state contained
% within "physics"

    figure(42)
		clf 

	%Hydrogen ion concentration
    subplot(3,3,1)
        physics.PlotNodal("H",-1, "Electrolyte");
        title("H^+")
		colorbar
	%oxygen concentration
	subplot(3,3,2)
        physics.PlotNodal("O2",-1, "Electrolyte");
        title("O_2")
		colorbar

	%FeOH+ concentration
	subplot(3,3,3)
        physics.PlotNodal("FeOH",-1, "Electrolyte");
        title("FeOH^{+}")
		colorbar
	
	%Electrolyte potential
	subplot(3,3,5)
        physics.PlotNodal("Epot",-1, "Electrolyte");
        title("$\varphi$",'Interpreter','latex')
		colorbar

	%Fe2+ concentration
	subplot(3,3,6)
        physics.PlotNodal("Fe",-1, "Electrolyte");
        title("Fe^{2+}")
		colorbar

	%Na+ concentration
	subplot(3,3,7)
        physics.PlotNodal("Na",-1, "Electrolyte");
        title("Na^+")
		colorbar

	%Cl- concentration
	subplot(3,3,8)
        physics.PlotNodal("Cl",-1, "Electrolyte");
        title("Cl^-")
		colorbar

	%OH- concentration
	subplot(3,3,9)
        physics.PlotNodal("OH",-1, "Electrolyte");
        title("OH^-")
		colorbar
	
	%total reaction currents and metal potential
	figure(45)
		clf
		yyaxis left
		plot(tvec/3600, I_an_vec)
		hold on
		plot(tvec/3600, I_cat_H_vec)
		plot(tvec/3600, I_cat_O_vec)
		xlabel('t [hours]')
		ylabel('$I_{anode} \;[A]$','Interpreter','latex')
		yyaxis right
		plot(tvec/3600, Em_vec)
		xlabel('t [hours]')
		ylabel('$E_m [V_{SHE}]$','Interpreter','latex')	
		legend('Corrosion','Hydrogen','Oxygen','Metal Potential')

	%Surface reaction rates
 	figure(44)
 		clf
 		physics.models{2}.plotReactions(physics);

	%other electrolyte-specific parameters (pH etc.)
	figure(43)
		clf
		physics.models{1}.plotFields(physics);

     drawnow();
end

