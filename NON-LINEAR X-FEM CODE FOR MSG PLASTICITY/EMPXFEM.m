% Main function of the 2D (plane strain) X-FEM code
% This version has the current capabilities: 3 material models, linear
% elements with 4 gauss points, displacement controlled loading.

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [SIGMA,S22,F1]=EMPXFEM(crack,XY,LE,MAT,PROP,Dnodes,SOL,...
    enrType,elem,graph,blend)

[NNodes, Ndof]=size(XY);
[NE, NNE] = size(LE);
[nfix, ~] = size(Dnodes);

% Identify and classify enriched elements and nodes
[typeEl,ElCrk,tipEl,splitEl,vertexEl,xTip,xVertex,enrNode]=...
    crackDetect(crack,XY,LE);

if strcmp(enrType,'geom')
 %geometrix/fixed area enrichment
 cnt=0;
 Renr=5e-4;
 xcent=crack(2,:);
 if graph==1 
  figure(1)
  theta = -pi:0.1:pi ;
  xo = xcent(1) + Renr*cos(theta) ;
  yo = xcent(2) + Renr*sin(theta) ;
  plot(xo,yo,'k-')
 end
 for in = 1:NNodes
  x=XY(in,:);
  dist=sqrt((x(1)-xcent(1))^2+(x(2)-xcent(2))^2)-Renr;
  if dist<0
   cnt=cnt+1;
   enrNode(in,1)=1;
  end
 end
end
% Initialize stiffness matrix and force vector
split=size(find(enrNode(:,1)==2),1);
tip=size(find(enrNode(:,1)==1),1);

% Ramp function definition
phiN=zeros(NNodes,1);
for in=1:NNodes
 if enrNode(in,1)>0
  phiN(in,1)=1;
 end
end

pos=zeros(NNodes,1);
nsnode=0;
ntnode=0;

for i=1:NNodes
 if (enrNode(i,1) == 2)
  pos(i,1)=(NNodes+nsnode*1+ntnode*4)+1;
  nsnode=nsnode+1;
 elseif (enrNode(i,1) == 1)  
  pos(i,1) = (NNodes + nsnode*1 + ntnode*4) + 1 ;
  ntnode = ntnode + 1 ;
 end
end
    
TU=NNodes*Ndof+split*1*Ndof+tip*4*Ndof;
    
if graph==1 
% Plot enriched nodes    
 for kj = 1:size(crack,1)-1
  cr = plot(crack(kj:kj+1,1),crack(kj:kj+1,2),'r-') ;
  set(cr,'LineWidth',3);
 end
 for kj = 1:size(crack,1)
  plot(crack(kj,1),crack(kj,2),'ro','MarkerFaceColor',[.49 1 .63],...
      'MarkerSize',5);
 end
  split_nodes = find(enrNode(:,1) == 2);
  tip_nodes   = find(enrNode(:,1) == 1);
  n1 = plot(XY(split_nodes,1),XY(split_nodes,2),'r*');
  n2 = plot(XY(tip_nodes,1),XY(tip_nodes,2),'rs');
  set(n1,'MarkerSize',15);
  set(n2,'MarkerSize',15);
end 

% Boundary conditions

isBdNode=false(TU,1);
Dirichlet(1:nfix)=Ndof*(Dnodes(1:nfix,1)-1) + Dnodes(1:nfix,2);
isBdNode(Dirichlet)=true;
bdNode=find(isBdNode);
freeNode=find(~isBdNode);
    
% Initialize the displacement and the incremental displacements
w = zeros(TU,1);
u=w;

% Initialize history dependent variables
SIGMA = zeros(3,3,300,NE);
EP=zeros(300,NE);
rG=zeros(300,NE);
EE=zeros(4,300,NE);
EPL=zeros(4,300,NE);
ETAP=zeros(300,NE);
stress_pnt=[];
strainp_val=[];
F1=[];
STRAINP=zeros(NE,4,4);

% Store solver controls
dt = 2./SOL(1); 

for step=1:SOL(1)

 loadfactor = step/SOL(1);
 err1 = 1.;    % Error/residual
 nit = 0;      % Number of iterations
 
 fprintf(1,'\n Step %f Load %f\n',step,loadfactor);
 
 while ((err1>SOL(2)) && (nit<SOL(3))) % Newton Raphson loop
   
  nit = nit + 1;
  % Compute the global stiffness matrix, force vector and residual 
  [GK,R,stress_pnt,STRAINP,strainp_val,F1]=globalKR(dt,Ndof,XY',NE,LE',...
    PROP,SIGMA,EP,w,MAT,EE,EPL,ETAP,NNE,TU,enrNode,ElCrk,typeEl,xTip,...
    xVertex,splitEl,tipEl,pos,vertexEl,step,nit,phiN,elem,blend,...
    stress_pnt,STRAINP,strainp_val,F1);

  % Prescribed displacements (elimination method)      
  dw=zeros(TU,1);
  dw(bdNode)=loadfactor*Dnodes(:,3)-w(bdNode)-u(bdNode);
  b=-R-GK*dw;

  % Solve system
  dw (freeNode)= GK(freeNode,freeNode)\b(freeNode);      

  % Check convergence
  w = w + dw;
  wnorm = dot(w,w);
  err1 = dot(dw,dw);
  err2 = dot(R(freeNode),R(freeNode));
  err1 = sqrt(err1/wnorm);
  err2 = sqrt(err2)/(Ndof*NNodes);
  fprintf(1,'Iteration number %d Correction %f Residual %f tolerance %f\n',nit,err1,err2,SOL(2));
  if nit == SOL(3)
   error('NEWTON-RAPHSON DID NOT CONVERGE');
  end    
 end
    
% Update the stress and accumulated plastic strain
 [SIGMA,EP,ETAP,EE,EPL,S22,STRAINP,strainp_val]=update_state(dt,Ndof,...
   XY',NE,LE',PROP,SIGMA,EP,w,MAT,ETAP,EE,EPL,NNE,enrNode,ElCrk,...
   typeEl,xTip,xVertex,splitEl,tipEl,vertexEl,pos,SOL,step,phiN,...
   graph,elem,blend,stress_pnt,STRAINP,strainp_val,F1);
     
% Update the total displacecment
 u=u+w;
end

% Emilio Martínez-Pañeda (mail@empaneda.com)