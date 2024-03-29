% Input file for the cracked plate example (Fig. 8, Comput Mech 2017)

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

% Initializing
close all
clear variables

% Flag variables
enrType='tip'; %Enrichment type: tip (topological) or geom (geometrical)
graph=0; %Graphical output: 1 activate, 0 deactive
blend=1; %Blending elements: 1 active, 0 deactive
elem='Q4'; %NOTE - THIS IS THE LINEAR ELEMENTS VERSION - ONLY WORKS FOR Q4

% Get nodal coordinates and element connectivity from ABAQUS mesh:
[XY,LE]=mesh2D;
% Crack definition
crack=[-14.1 0;0 0];

% Plot the mesh - including the crack
if graph==1
 plotmesh(XY,LE,'b-',elem)
 for i=1:size(crack,1)-1
  cr=plot(crack(i:i+1,1),crack(i:i+1,2),'r-');
  set(cr,'LineWidth',3);
 end
 for i=1:size(crack,1)
  plot(crack(i,1),crack(i,2),'ro','MarkerFaceColor',...
      [.49 1 .63],'MarkerSize',5);
 end
end

% Material model definition
% MAT=0 - Linear elasticity - PROP=[YOUNG NU]
% MAT=1 - J2 plasticity (without hardening) - PROP=[YOUNG NU SYIELD] 
% MAT=2 - CMSG plasticity - PROP=[YOUNG NU SYIELD L N]
MAT=2;
PROP=[260000 0.3 200 0.005 0.2];

% Prescribed displacements [Node, DOF, Magnitude]
[NUMNP, ~] = size(XY);
j=1;opt=-100;
for i=1:size(XY,1)
 if XY(i,2)<-49.999
  Dnodes(j,:)=[i 2 -0.0011]; % U top
  j=j+1;
 end
 if XY(i,2)>49.999
  Dnodes(j,:)=[i 2 0.0011]; % U bottom
  j=j+1;
 end    
 if XY(i,1)>20.999
  if XY(i,2)<0
   if XY(i,2)>opt
    opt=XY(i,2); 
    pos=i;
   end
  end
 end
end
Dnodes(j,:)=[pos 1 0]; % Constraint rigid body motion
Dnodes=sortrows(Dnodes,1);

% Solution controls [Num. of increments, Conv. tolerance, max. iterations]
SOL=[20, 1E-5, 50];

% Calling the main function and storing the stress tensor
[SIGMA,S22,F1]=EMPXFEM(crack,XY,LE,MAT,PROP,Dnodes,SOL,enrType,...
    elem,graph,blend);

% Plotting results (opening stress distribution ahead of the crack)
coord=[0.00006 0; 0.00007 0; 0.00008 0; 0.00009 0; 0.00015 0; 0.0001 0; 
  0.0002 0; 0.0003 0; 0.0004 0; 0.0005 0; 0.0006 0; 0.0007 0; 0.0008 0;
  0.0009 0; 0.001 0; 0.002 0; 0.003 0; 0.004 0; 0.005 0; 0.006 0; 
  0.007 0; 0.008 0; 0.009 0; 0.01 0; 0.015 0; 0.02 0; 0.025 0; 0.03 0;
  0.035 0; 0.04 0; 0.05 0; 0.06 0; 0.07 0; 0.08 0; 0.09 0; 0.1 0; 0.12 0;
  0.14 0; 0.16 0; 0.18 0; 0.2 0; 0.24 0; 0.28 0; 0.32 0; 0.38 0; 0.44 0; 
  0.5 0; 0.6 0; 0.7 0; 0.8 0; 0.9 0; 1 0; 1.5 0; 2 0; 3 0;5 0; 10 0; 19 0];
F1.Values=S22;
vq=F1(coord(:,1),coord(:,2));
fileID = fopen('FEMref.txt','r');
formatSpec = '%f %f';
sizeA = [2 Inf];
Ab = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
figure
x=Ab(1,:)/PROP(4);
y=Ab(2,:)/PROP(3);
p=semilogx(x,y,'-k');
p(1).LineWidth = 2;
xlabel('$r/l$','Interpreter','LaTex','FontSize',20) 
ylabel('$\sigma_{22}/\sigma_{Y}$','Interpreter','LaTex','FontSize',20)
xlim([0.002 20])
ylim([0 15])
xlim manual
hold on
p1=semilogx(coord(:,1)/PROP(4),vq/PROP(3),'gs','MarkerSize',10,'MarkerFaceColor','g');
legend('FEM','X-FEM')

% Emilio Martínez-Pañeda (mail@empaneda.com)
