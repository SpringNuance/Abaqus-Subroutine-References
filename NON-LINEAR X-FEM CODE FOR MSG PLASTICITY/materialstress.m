% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function stress = materialstress(ndof,strain,materialprops,stressN)
STRESS(1)=stressN(1);
STRESS(2)=stressN(2);
STRESS(3)=stressN(4);
STRESS(4)=stressN(3);
C = ELAST(materialprops);
stress = zeros(ndof,ndof);
stressT=C*strain;
stress(1,1) = stressT(1)+STRESS(1);
stress(2,2) = stressT(2)+STRESS(2);
stress(1,2) = stressT(3)+STRESS(4);
stress(2,1) = stress(1,2);