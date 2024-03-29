close all; clear all; clc

fileName = "Results/Em_0/10.mat";
fileTitle = "$E_m = 0\;\mathrm{V}_{SHE}$"

addpath(genpath('./Models'))
addpath(genpath('./Shapes'))

load(fileName,'tvec','CL_vec','Cmax_vec','mesh','physics')	

%% Time vs avarage and maximum hydrogen contents
f1 = figure;
	yyaxis left

	% Average concentrations are saved within CL_vec
	plot(tvec/3600/24,CL_vec,'LineWidth',2, 'DisplayName', fileTitle);
	hold on
	xlabel('$time \;[\mathrm{days}]$','Interpreter','latex')
	ylabel('$\overline{C_L} \;[\mathrm{mol}/\mathrm{m}^3]$','Interpreter','latex')
	yyaxis right

	%maximum concentrations saved within Cmax_vec
	plot(tvec/3600/24,Cmax_vec,'LineWidth',2, 'DisplayName', fileTitle);
	ylabel('$\tilde{C_L} \;[\mathrm{mol}/\mathrm{m}^3]$','Interpreter','latex')

	savefig(f1, "Figures/HydrogenOverTime")

%% Surface plot of pH and interstitial lattice hydrogen contents
f2 = figure;
t = tiledlayout(1,1);

	%Interstitial lattice hydrogen, plotted at nodal values
	ax1 = axes(t);
	physics.PlotNodal("CL",-1, "Metal");
	xlim([-0.01 0.01])
	ylim([0 0.01])
	colormap(ax1,'cool')
	cb1 = colorbar;
	axis off 

	%pH, plotted through a post-processing function contained within the
	%electrolyte model
	ax2 = axes(t);
	physics.models{7}.plotpH(physics);
	colormap(ax2,'jet')
	cb2 = colorbar;
	xlim([-0.01 0.01])
	ylim([0 0.01])
	axis off 

	%combining the two plots
	cb1.Title.String = {'$C_L$', '[$\mathrm{mol}/\mathrm{m}^3$]', ' '};
	cb1.Title.Interpreter='latex';
	
	cb2.Title.String = {'pH', '[$-$]', ' '};
	cb2.Title.Interpreter='latex';
	
	fg = gcf;
	fg = gcf;
	fg.Units = 'centimeters';
	fg.Position(3) = 12;
	fg.Position(4) = 7;
	
	t.InnerPosition = [0.1,0.1,0.8,0.8];
	t.OuterPosition = [0.02,0,0.94,1];
	cb1.Position = [0.91 0.1318 0.0306 0.6];
	cb2.Position = [0.05 0.1318 0.0306 0.6];
	drawnow()
	fg.Position(4) = 6.5;
	drawnow()

	savefig(f2, "Figures/SurfacePlot")

%% Surface reaction rates
f3 = figure;
	physics.models{9}.plotReactions(physics);
	savefig(f3, "Figures/ReactionRates")

%% Surface reaction rates
f4 = figure;
	physics.PlotNodal("Theta",-1, "Interface");
	cb = colorbar;
	cb.Title.String = {'$\theta$', '[$\mathrm{mol}/\mathrm{m}^2$]', ' '};
	cb.Title.Interpreter='latex';
	savefig(f4, "Figures/SurfaceCoverage")

function savefig(fg, savepath)
	print(fg, savepath+".png",'-dpng','-r1200')
	print(fg, savepath+".jpg",'-djpeg','-r1200')
    print(fg, savepath+".eps",'-depsc','-r1200')
end
