% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [stress1,ep1,r1,etap1,ee1,epl1,stress_val,strainp_val,Icount]=...
   update_el_state(dt,ndof,materialprops,stress0,eplas0,U,MAT,nelem,...
   egrad0,ee0el,epl0el,type_elem,enrich_node,elem_crk,xVertex,node,...
   element,W,Q,tip_elem,SOL,step,stress_val,coordN,phiN,elem,blend,...
   STRAINP,strainp_val,Icount)

if MAT>0   
 strain=zeros(4,1);
else
 strain=zeros(3,1);
end
stressi0=zeros(4,1);
ee0=zeros(4,1);
epl0=zeros(4,1);
       
for kk=1:size(W,1)
 B=[];
 Gpt=Q(kk,:);
 B=[B xfemBmat(Gpt,nelem,type_elem,enrich_node(:,1),elem_crk,xVertex,...
     1,node,element,MAT,tip_elem,phiN,elem,blend)] ;
 strain=B*U;       
    
 if MAT>0
  strain(4)=0.0; % Plane strain
  ep0=eplas0(kk);
  eta0=egrad0(kk);
 end
 
 for j=1:4
  stressi0(j)=stress0(j,kk);
  ee0(j)=ee0el(j,kk);
  epl0(j)=epl0el(j,kk);
 end

 if MAT==2 
  for j=1:3
   deps(j)=strain(j);
  end 
  [stress,dep,dee,depl,detap,strainp_val,Icount]=...
      combHardC(materialprops,deps',ep0,ee0,epl0,eta0,nelem,coordN,...
      Gpt(1),Gpt(2),STRAINP,strainp_val,Icount);
 elseif MAT==1
  dep=deplas(dt,stressi0,ep0,strain,materialprops);
  stress=materialstressP(stressi0,dep,strain,materialprops);
 else
  stress=materialstress(ndof,strain,materialprops,stressi0);
 end
      
 if step==SOL(1) %Plot stress contours
  stress_val = [stress_val; stress(2,2)] ;      
 end    
      
 if MAT > 0   
  for i=1:3
   for j=1:3
    stress1(i,j,kk)=stress(i,j);
   end
  end
 else
  ep0=0;
  dep=0;
  dr=0;
  dee=0;
  depl=0;
  detap=0;
  for i=1:2
   for j=1:2
    stress1(i,j,kk) = stress(i,j);
   end
  end
 end

 ep1(kk) = ep0 + dep;
 if MAT == 2
  dr=0;
  etap1(kk)=detap;
  ee1(:,kk)=dee;
  epl1(:,kk)=depl;
 else
  etap1=zeros(kk); ee1=zeros(4,kk); epl1=zeros(4,kk);
 end
 r1(kk)= dr;     
end