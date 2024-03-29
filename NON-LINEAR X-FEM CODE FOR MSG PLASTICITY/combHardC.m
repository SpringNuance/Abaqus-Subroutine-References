% Function for CMSG plasticity

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function[stress,dp,EELAS,EPLAS,EQGRAD,strainp_val,Icount]=combHardC(mp,...
    deps,epN,eeN,eplN,etaN,noel,coordN,s,t,STRAINP,strainp_val,Icount)

mu=mp(1)/(2*(1+mp(2)));
K=mp(1)/(3*(1-2*mp(2)));
if any(deps)==0
 stress=zeros(3,3); 
 dp=0.0; 
 EELAS=zeros(1,4); 
 EPLAS=zeros(1,4); 
 EQGRAD=0.0;
 Icount=1;
else 
 Iden=[1 1 0 1]'; 
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
 sigmae=mp(3);
 h=0;rhs=mp(3);
 while abs(rhs)>1E-6*mp(3)
  rhs=3*mu*(def-defi*((sigmae/SigmaF)^20))-sigmae;
  sigmae=sigmae+rhs/(3*mu*h+1);
  h=20*defi*((sigmae/SigmaF)^19)*(1/SigmaF);
 end
 dp=def-sigmae/(3*mu);
 dstr=strain*2*mu/(1+3*mu*dp/sigmae);
 dpstrn=dstr*(3*dp)/(2*sigmae);
 dpstrn(3)=2*dpstrn(3);
 destrn=deps-dpstrn; 
 EELAS=eeN+destrn;
 EPLAS=eplN+dpstrn;
 I3=EELAS(1)+EELAS(2)+EELAS(4);
 stress(1,1)=dstr(1)+K*I3;
 stress(2,2)=dstr(2)+K*I3;
 stress(1,2)=dstr(3);
 stress(3,3)=dstr(4)+K*I3;
 stress(2,1)=stress(1,2);
 stress(1,3)=0.0;
 stress(2,3)=0.0;
 stress(3,1)=stress(1,3);
 stress(3,2)=stress(2,3);
 
 strainp_val(Icount,:)=dpstrn';
 Icount=Icount+1;
 
 deriv(1,1)=-1/4*(1-t);
 deriv(1,2)=1/4*(1-t);
 deriv(1,3)=1/4*(1+t);
 deriv(1,4)=-1/4*(1+t);
 deriv(2,1)=-1/4*(1-s); 
 deriv(2,2)=-1/4*(1+s);
 deriv(2,3)=1/4*(1+s);
 deriv(2,4)=1/4*(1-s);
 
 xjacm(1,1)=deriv(1,1)*coordN(1,1)+deriv(1,2)*coordN(2,1)+deriv(1,3)*...
     coordN(3,1)+deriv(1,4)*coordN(4,1);
 xjacm(1,2)=deriv(1,1)*coordN(1,2)+deriv(1,2)*coordN(2,2)+deriv(1,3)*...
     coordN(3,2)+deriv(1,4)*coordN(4,2);
 xjacm(2,1)=deriv(2,1)*coordN(1,1)+deriv(2,2)*coordN(2,1)+deriv(2,3)*...
     coordN(3,1)+deriv(2,4)*coordN(4,1);
 xjacm(2,2)=deriv(2,1)*coordN(1,2)+deriv(2,2)*coordN(2,2)+deriv(2,3)*...
     coordN(3,2)+deriv(2,4)*coordN(4,2);

 djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1);
 xjaci(1,1)=xjacm(2,2)/djacb;
 xjaci(1,2)=-xjacm(1,2)/djacb;
 xjaci(2,1)=-xjacm(2,1)/djacb;
 xjaci(2,2)=xjacm(1,1)/djacb;
 
 A1=xjaci(1,1)*deriv(1,1)+xjaci(1,2)*deriv(2,1);
 A2=xjaci(1,1)*deriv(1,2)+xjaci(1,2)*deriv(2,2);
 A3=xjaci(1,1)*deriv(1,3)+xjaci(1,2)*deriv(2,3);
 A4=xjaci(1,1)*deriv(1,4)+xjaci(1,2)*deriv(2,4);
 B1=xjaci(2,1)*deriv(1,1)+xjaci(2,2)*deriv(2,1);
 B2=xjaci(2,1)*deriv(1,2)+xjaci(2,2)*deriv(2,2);
 B3=xjaci(2,1)*deriv(1,3)+xjaci(2,2)*deriv(2,3);
 B4=xjaci(2,1)*deriv(1,4)+xjaci(2,2)*deriv(2,4);
 
 eta(1)=A1*STRAINP(noel,1,1)+A2*STRAINP(noel,2,1)+A3*STRAINP(noel,3,1)...
     +A4*STRAINP(noel,4,1);
 eta(2)=A1*STRAINP(noel,1,3)+A2*STRAINP(noel,2,3)+A3*STRAINP(noel,3,3)...
     +A4*STRAINP(noel,4,3)-B1*STRAINP(noel,1,1)-B2*STRAINP(noel,2,1)...
     -B3*STRAINP(noel,3,1)-B4*STRAINP(noel,4,1);
 eta(3)=B1*STRAINP(noel,1,1)+B2*STRAINP(noel,2,1)+B3*STRAINP(noel,3,1)...
     +B4*STRAINP(noel,4,1);
 eta(4)=A1*STRAINP(noel,1,2)+A2*STRAINP(noel,2,2)+A3*STRAINP(noel,3,2)...
     +A4*STRAINP(noel,4,2);
 eta(5)=A1*STRAINP(noel,1,4)+A2*STRAINP(noel,2,4)+A3*STRAINP(noel,3,4)...
     +A4*STRAINP(noel,4,4);
 eta(6)=B1*STRAINP(noel,1,1)+B2*STRAINP(noel,2,1)+B3*STRAINP(noel,3,1)...
     +B4*STRAINP(noel,4,1);
 eta(7)=A1*STRAINP(noel,1,2)+A2*STRAINP(noel,2,2)+A3*STRAINP(noel,3,2)...
     +A4*STRAINP(noel,4,2);
 eta(8)=B1*STRAINP(noel,1,3)+B2*STRAINP(noel,2,3)+B3*STRAINP(noel,3,3)...
     +B4*STRAINP(noel,4,3)-A1*STRAINP(noel,1,2)-A2*STRAINP(noel,2,2)...
     -A3*STRAINP(noel,3,2)-A4*STRAINP(noel,4,2);
 eta(9)=B1*STRAINP(noel,1,2)+B2*STRAINP(noel,2,2)+B3*STRAINP(noel,3,2)...
     +B4*STRAINP(noel,4,2);
 eta(10)=B1*STRAINP(noel,1,4)+B2*STRAINP(noel,2,4)+B3*STRAINP(noel,3,4)...
     +B4*STRAINP(noel,4,4); 
 eta(11)=A1*STRAINP(noel,1,4)+A2*STRAINP(noel,2,4)+A3*STRAINP(noel,3,4)...
     +A4*STRAINP(noel,4,4);
 eta(12)=B1*STRAINP(noel,1,4)+B2*STRAINP(noel,2,4)+B3*STRAINP(noel,3,4)...
     +B4*STRAINP(noel,4,4);
 eta(13)=-(A1*STRAINP(noel,1,4)+A2*STRAINP(noel,2,4)...
     +A3*STRAINP(noel,3,4)+A4*STRAINP(noel,4,4));
 eta(14)=-(B1*STRAINP(noel,1,4)+B2*STRAINP(noel,2,4)...
     +B3*STRAINP(noel,3,4)+B4*STRAINP(noel,4,4));
 
 etat=eta(1)^2+eta(2)^2+eta(3)^2+eta(4)^2+eta(5)^2+eta(6)^2+eta(7)^2 ...
     +eta(8)^2+eta(9)^2+eta(10)^2+eta(11)^2+eta(12)^2+eta(13)^2+eta(14)^2;
 deqgrad=sqrt(1/4*etat);
 EQGRAD=etaN+deqgrad; 
 end 
end