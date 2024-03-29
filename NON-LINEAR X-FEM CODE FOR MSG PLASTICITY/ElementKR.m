% Function to compute the element stiffness matrix

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [dsde,stresst,strainp_val,Icount]=ElementKR(dt,ndof,...
    materialprops,stress0,eplas0,MAT,egrad0,ee0el,epl0el,nelem,intpt,...
    strain,coordN,xi,STRAINP,strainp_val,Icount)

stressi0=zeros(4,1);
ee0=zeros(4,1);
epl0=zeros(4,1);
      
if MAT > 0
 strain(4)=0.0; % Plane strain
 ep0=eplas0(intpt);
end

for j=1:4
 stressi0(j)=stress0(j,intpt);
 ee0(j)=ee0el(j,intpt);
 epl0(j)=epl0el(j,intpt);
end

if MAT == 2
 for j=1:3
  deps(j) = strain(j);
 end  
 ep0=eplas0(intpt);
 eta0=egrad0(intpt);
 D=ELAST(materialprops);
 dsde=combHardTanC(materialprops,D,deps',ep0,ee0,eta0);
 [stress,~,~,~,~,strainp_val,Icount]=combHardC(materialprops,deps',ep0,...
     ee0,epl0,eta0,nelem,coordN,xi(1),xi(2),STRAINP,strainp_val,Icount);
elseif MAT==1
 ep0=eplas0(intpt);
 dep=deplas(dt,stressi0,ep0,strain,materialprops);
 stress=materialstressP(stressi0,dep,strain,materialprops);
 dsde=materialstiffness(stress,ep0,dep,strain,materialprops);
else
 dsde=ELAST(materialprops);
 stress=materialstress(ndof,strain,materialprops,stressi0);	
end
stresst(1)=stress(1,1);
stresst(2)=stress(2,2);
stresst(3)=stress(1,2);   
