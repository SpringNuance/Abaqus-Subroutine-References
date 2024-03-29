% Function for CMSG plasticity

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [Dtan]=combHardTanC(mp,D,deps,epN,eeN,etaN)
mu=mp(1)/(2*(1+mp(2)));
K=mp(1)/(3*(1-2*mp(2)));
if any(deps)==0
 Dtan = D;
else
 Iden = [1 1 0 1]'; 
 I1=eeN(1)+eeN(2)+eeN(4);
 strad=eeN-I1*Iden/3;
 deps(4)=0;
 I2=deps(1)+deps(2)+deps(4);
 dstrad=deps-I2*Iden/3;
 strain=dstrad+strad;
 strain(3)=strain(3)/2;
 dstrad(3)=dstrad(3)/2;
 def=sqrt((2/3)*(strain(1)^2+strain(2)^2+strain(4)^2+2*(strain(3)^2)));
 defi=sqrt((2/3)*(dstrad(1)^2+dstrad(2)^2+dstrad(4)^2+2*(dstrad(3)^2)));
 SigmaF=mp(3)*((mp(1)/mp(3))^mp(5))*sqrt((epN+mp(3)/mp(1))^(2*mp(5))...
     +mp(4)*etaN);
 sigmae=mp(3);h=0;rhs=mp(3);
 while abs(rhs)>1E-6*mp(3)
  rhs=3*mu*(def-defi*((sigmae/SigmaF)^20))-sigmae;
  sigmae=sigmae+rhs/(3*mu*h+1);
  h=20*defi*((sigmae/SigmaF)^19)*(1/SigmaF);
 end
 dp=def-sigmae/(3*mu);
 dstr=strain*2*mu/(1+3*mu*dp/sigmae);
 Q=2/3*sigmae/def;
 R=((h-dp/sigmae)/(def*sigmae))*(3*mu/(1+3*mu*h));
 Dtan(1,1)=Q+(K-Q*1/3)-R*dstr(1)*dstr(1);
 Dtan(2,2)=Q+(K-Q*1/3)-R*dstr(2)*dstr(2);
 Dtan(1,2)=(K-Q*1/3)-R*dstr(1)*dstr(2);
 Dtan(2,1)=Dtan(1,2);
 Dtan(2,3)=-R*dstr(2)*dstr(3);
 Dtan(3,2)=Dtan(2,3); 
 Dtan(1,3)=-R*dstr(1)*dstr(3);
 Dtan(3,1)=Dtan(1,3);
 Dtan(3,3)=Q/2-R*dstr(3)*dstr(3);            
 end
end