% Function to compute the branch functions
% Based on the original function by Nguyen Vinh Phu

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [Br,dBdx,dBdy]=branch_gp(r,theta,alpha,MAT)
if( r ~=0 )
 if MAT==2 
  r2=r^0.36163; % see Shi et al. (2001)
 else
  r2=sqrt(r);
 end
else
 r2=0.1d-4;
 theta=0.0d0 ;
end

if MAT==2
 fac=(0.36163)*(r^(-0.63837)); % see Shi et al. (2001)
else
 fac=0.5/r2 ;
end

st2=sin(theta/2.);
ct2=cos(theta/2.);
s3t2=sin(1.5*theta);
c3t2=cos(1.5*theta);
st=sin(theta);
ct=cos(theta);

% LEFM Functions 
Br(1)=r2*st2;
Br(2)=r2*ct2;
Br(3)=r2*st2*st;
Br(4)=r2*ct2*st;

% Derivatives in local coords (x1,x2)
dPhi1dx1=-fac*st2;
dPhi1dx2=fac*ct2;
dPhi2dx1=dPhi1dx2;
dPhi2dx2=-dPhi1dx1;
dPhi3dx1=-fac*s3t2*st;
dPhi3dx2=fac*(st2+s3t2*ct);
dPhi4dx1=-fac*c3t2*st;
dPhi4dx2=fac*(ct2+c3t2*ct);

% Derivatives in global coords (x,y)

dx1dx=cos(alpha); 
dx2dx = -sin(alpha);
dx1dy=sin(alpha); 
dx2dy = cos(alpha);

dBdx(1)=dPhi1dx1*dx1dx+dPhi1dx2*dx2dx;
dBdy(1)=dPhi1dx1*dx1dy+dPhi1dx2*dx2dy;
dBdx(2)=dPhi2dx1*dx1dx+dPhi2dx2*dx2dx;
dBdy(2)=dPhi2dx1*dx1dy+dPhi2dx2*dx2dy;
dBdx(3)=dPhi3dx1*dx1dx+dPhi3dx2*dx2dx;
dBdy(3)=dPhi3dx1*dx1dy+dPhi3dx2*dx2dy;
dBdx(4)=dPhi4dx1*dx1dx+dPhi4dx2*dx2dx;
dBdy(4)=dPhi4dx1*dx1dy+dPhi4dx2*dx2dy;
