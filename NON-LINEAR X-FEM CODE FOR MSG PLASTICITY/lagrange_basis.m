% Function to compute the lagrange interpolant basis
% Based on the original function by Jack Chessa

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [Nv,dNdxi]=lagrange_basis(type,coord,flag)
 
if(nargin==2)
 flag=1;
end
switch type  
 case 'T3'
  if size(coord,2) < 2
    disp('Error two coordinates needed for the T3 element')
  else
   xi=coord(1); eta=coord(2);
   N=[1-xi-eta;xi;eta];
   dNdxi=[-1,-1;1,0;0,1];
  end
    
  case 'T3fs'
   if size(coord,2) < 2
    disp('Error two coordinates needed for the T3fs element')
   else
    xi=coord(1); eta=coord(2);
    N=[1-xi-eta;xi;eta];
    dNdxi=[-1,-1;1,0;0,1];
   end    
 case 'Q4'
  if size(coord,2)<2
   disp('Error two coordinates needed for the Q4 element')
  else
   xi=coord(1); 
   eta=coord(2);
   N=1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
   dNdxi=1/4*[-(1-eta),-(1-xi);1-eta,-(1+xi);1+eta,1+xi;-(1+eta),1-xi];
  end
 case 'Q8'
  if size(coord,2)<2
   disp('Error two coordinates needed for the Q8 element')
  else
   xi=coord(1); 
   eta=coord(2);
   N=[-0.25*(1-xi)*(1-eta)*(1+xi+eta);-0.25*(1+xi)*(1-eta)*(1-xi+eta);...
     -0.25*(1+xi)*(1+eta)*(1-xi-eta);-0.25*(1-xi)*(1+eta)*(1+xi-eta);...
     0.5*(1-xi^2)*(1-eta);0.5*(1+xi)*(1-eta^2);0.5*(1-xi^2)*(1+eta);...
     0.5*(1-xi)*(1-eta^2)];
   dNdxi=[0.25*(1-eta)*(2*xi+eta), 0.25*(1-xi)*(xi+2*eta);...
        0.25*(1-eta)*(2*xi-eta), -0.25*(1+xi)*(xi-2*eta);...
        0.25*(1+eta)*(2*xi+eta), 0.25*(1+xi)*(xi+2*eta);...
        0.25*(1+eta)*(2*xi-eta), -0.25*(1-xi)*(xi-2*eta);...
        -(1-eta)*xi, -0.5*(1-xi^2);0.5*(1-eta^2), -(1+xi)*eta;...
        -(1+eta)*xi, 0.5*(1-xi^2);-0.5*(1-eta^2), -(1-xi)*eta]; 
   end    
  otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
 end
 
  I=eye(flag);
  Nv=[];
  for i=1:size(N,1)
   Nv=[Nv;I*N(i)];
  end
  
if (flag==1)
 B=dNdxi;
elseif(flag==2)
 B=zeros(flag*size(N,1),3);
 B(1:flag:flag*size(N,1)-1,1)=dNdxi(:,1);
 B(2:flag:flag*size(N,1),2)=dNdxi(:,2);
 B(1:flag:flag*size(N,1)-1,3)=dNdxi(:,2);
 B(2:flag:flag*size(N,1),3)=dNdxi(:,1);
elseif (flag==3)
 B=zeros(flag*size(N,1),6);   
 disp('Error: need to add 3D N and dNdxi')
 B(1:flag:flag*size(N,1)-2,1)=dNdxi(:,1);
 B(2:flag:flag*size(N,1)-1,2)=dNdxi(:,2);
 B(3:flag:flag*size(N,1),3)=dNdxi(:,3);
 B(2:flag:flag*size(N,1)-1,4)=dNdxi(:,3);
 B(3:flag:flag*size(N,1),4)=dNdxi(:,2);
 B(3:flag:flag*size(N,1),5)=dNdxi(:,1);
 B(1:flag:flag*size(N,1)-2,5)=dNdxi(:,3);
 B(1:flag:flag*size(N,1)-2,6)=dNdxi(:,2);
 B(2:flag:flag*size(N,1)-1,6)=dNdxi(:,1);
end 