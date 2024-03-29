% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function S=materialstressP(stress,dep,dstrain,mp)

K=mp(1)/(3*(1-2*mp(2)));
E=mp(1);
nu=mp(2);

S=zeros(3,3);
de=zeros(3,3);
dl=[[1,0,0];[0,1,0];[0,0,1]];  

devol=trace(dstrain);
p=trace(stress);
   
se=0.;
for i=1:3
 for j=1:3
  de(i,j)=dstrain(i,j)-dl(i,j)*devol/3;
  S(i,j)=stress(i,j)-dl(i,j)*p/3+E/(1+nu)*de(i,j);
  se=se+S(i,j)*S(i,j);
 end
end
se=sqrt(1.5*se);

if(se>0)
 beta=1-1.5*E*dep/((1+nu)*se);
else
 beta=1.;
end            
for i=1:3
 for j=1:3
  S(i,j)=beta*S(i,j)+(p/3+K*devol)*dl(i,j);
 end
end