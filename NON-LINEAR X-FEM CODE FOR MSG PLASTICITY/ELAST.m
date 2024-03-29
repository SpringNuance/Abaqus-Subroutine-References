% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function C=ELAST(PROP)
C = zeros(3,3);
C(1,1)=(1-PROP(2))*PROP(1)/((1+PROP(2))*(1-2*PROP(2)));
C(2,2)=(1-PROP(2))*PROP(1)/((1+PROP(2))*(1-2*PROP(2)));
C(1,2)=PROP(2)*PROP(1)/((1+PROP(2))*(1-2*PROP(2)));
C(2,1)=C(1,2);
C(1,3)=0;
C(2,3)=0;
C(3,1)=C(1,3);
C(2,3)=C(3,2);
C(3,3)=PROP(1)/(2*(1+PROP(2)));