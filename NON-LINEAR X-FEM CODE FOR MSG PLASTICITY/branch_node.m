% Function to compute the branch functions
% Based on the original function by Nguyen Vinh Phu

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [BrI]=branch_node(r,theta,MAT)

if(r ~=0)
 if MAT==2
  r2=r^0.36163; % see Shi et al. (2001)
 else
  r2=sqrt(r);
  end
else
 r2=0.1d-4;
 theta=0.0d0;
end
st2=sin(theta/2.);
ct2=cos(theta/2.);
st=sin(theta);

% Functions 
BrI(1)=r2*st2;
BrI(2)=r2*ct2;
BrI(3)=r2*st2*st;
BrI(4)=r2*ct2*st;
        