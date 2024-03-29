% Function to compute the plastic strain increment for a given strain inc.
% Based on the original function by Alan Bower

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function dep = deplas(dt,stress,eplas,dstrain,mp)

E=mp(1);
nu=mp(2);
Y=mp(3);
e0=mp(4);
n=mp(5);
edot0=mp(6);
m=mp(7);

S=zeros(3,3);
de=zeros(3,3);
dl=[[1,0,0];[0,1,0];[0,0,1]];  

devol=trace(dstrain);
p=trace(stress);
   
sequiv=0.;
for i=1:3
 for j=1:3
  de(i,j)=dstrain(i,j)-dl(i,j)*devol/3;
  S(i,j)=stress(i,j)-dl(i,j)*p/3+E/(1+nu)*de(i,j);
  sequiv=sequiv+S(i,j)*S(i,j);
 end
end
sequiv=sqrt(1.5*sequiv);

dep=10^(-15);
err=Y;
tol=10^(-06)*Y;
if(sequiv*edot0==0)
 dep = 0.;
else
 while (err>tol)
  c=(1+(eplas+dep)/e0)^(1/n)*(dep/(dt*edot0))^(1/m);
  f=sequiv/Y-1.5*dep*E/(Y*(1+nu))-c;
  dfde=-1.5*E/(Y*(1+nu))-c*(1/(n*(eplas+dep+e0))+1/(m*dep));
  enew=dep-f/dfde;
  if (enew<0)
   dep=dep/10;
  else
   dep=enew;
  end
   err=abs(f);
 end
end
