% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [stress_new,eplas_new,eta_new,ee_new,ep_new,stress_val,STRAINP,...
 strainp_val]=update_state(dt,ndof,coords,nelem,connect,materialprops,...
 stress,eplas,dofs,MAT,eta,ee,ep,nne,enrich_node,elem_crk,type_elem,...
 xTip,xVertex,split_elem,tip_elem,vertex_elem,pos,SOL,step,phiN,...
 graph,elem,blend,stress_pnt,STRAINP,strainp_val,F1)

node=coords';
element=connect';
[NCorN, ~]=size(node);

lmndof=zeros(ndof,nne);
   
stress_val=[];
STRAINPT=zeros(NCorN,4);
Icount=1;
   
if MAT==2
 F1.Values=strainp_val(:,1);
 STRAINPT(:,1) = F1(node);
 F1.Values=strainp_val(:,2); 
 STRAINPT(:,2) = F1(node);
 F1.Values=strainp_val(:,3);
 STRAINPT(:,3) = F1(node);
 F1.Values=strainp_val(:,4);
 STRAINPT(:,4) = F1(node); 
 strainp_val=zeros(size(stress_pnt,1),4);
end  

for lmn=1:nelem
 sctr=element(lmn,:);
 coordN=node(sctr',:);
 
 [Nlines, ~]=size(coordN);
 if Nlines>4
  coordN(5,:)=[];  coordN(5,:)=[];  coordN(5,:)=[]; coordN(5,:)=[];
 end 
 
 if MAT==2
  sctrT=sctr(1,1:4);
  STRAINP(lmn,:,1)=STRAINPT(sctrT',1);
  STRAINP(lmn,:,2)=STRAINPT(sctrT',2);
  STRAINP(lmn,:,3)=STRAINPT(sctrT',3);  
  STRAINP(lmn,:,4)=STRAINPT(sctrT',4);     
 end

 [W,Q]=gauss_rule(lmn,enrich_node,elem_crk,xTip,xVertex,tip_elem,...
    split_elem,vertex_elem,node,element,elem);

 U=[];
 U=[U; element_disp(lmn,pos(:,1),enrich_node(:,1),dofs,element)];

 for a=1:nne
  for i=1:ndof
   lmndof(i,a)=dofs(ndof*(connect(a,lmn)-1)+i);
  end
 end
 nintp=size(W,1);
 lmnstress=zeros(4,nintp);
 lmneplas=zeros(nintp,1);
 lmnEPL=zeros(4,nintp);
 lmnEE=zeros(4,nintp);
 lmnETAP=zeros(nintp,1);      
 if MAT>0   
  for a=1:nintp
   lmnstress(1,a)=stress(1,1,a,lmn);
   lmnstress(2,a)=stress(2,2,a,lmn);
   lmnstress(4,a)=stress(3,3,a,lmn);
   lmnstress(3,a)=stress(1,2,a,lmn);
   lmneplas(a)=eplas(a,lmn);
   lmnEPL(1,a)=ep(1,a,lmn);
   lmnEPL(2,a)=ep(2,a,lmn);
   lmnEPL(3,a)=ep(3,a,lmn);
   lmnEPL(4,a)=ep(4,a,lmn);
   lmnEE(1,a)=ee(1,a,lmn);
   lmnEE(2,a)=ee(2,a,lmn);
   lmnEE(3,a)=ee(3,a,lmn);
   lmnEE(4,a)=ee(4,a,lmn);
   lmnETAP(a)=eta(a,lmn);   
  end
 end
      
 [lmnstress,lmneplas,lmnR,lmnETAP,lmnEE,lmnEPL,stress_val,strainp_val,...
     Icount]=update_el_state(dt,ndof,materialprops,lmnstress,lmneplas,U,...
     MAT,lmn,lmnETAP,lmnEE,lmnEPL,type_elem,enrich_node,elem_crk,...
     xVertex,node,element,W,Q,tip_elem,SOL,step,stress_val,coordN,phiN,...
     elem,blend,STRAINP,strainp_val,Icount);
 
 if MAT>0   
  for a=1:nintp
   for i=1:3
    for j=1:3
     stress_new(i,j,a,lmn)=lmnstress(i,j,a);
    end
   end
   for j=1:4
    ee_new(j,a,lmn)=lmnEE(j,a);
    ep_new(j,a,lmn)=lmnEPL(j,a);
   end
   eplas_new(a,lmn)=lmneplas(a);
   r_new(a,lmn)=lmnR(a);
   eta_new(a,lmn)=lmnETAP(a);
  end    
 else
  for a=1:nintp
   for i=1:2
    for j=1:2
     stress_new(i,j,a,lmn)=lmnstress(i,j,a);
    end
   end
   for j=1:4
    ee_new(j,a,lmn)=lmnEE(j,a);
    ep_new(j,a,lmn)=lmnEPL(j,a);
   end        
   eplas_new(a,lmn)=lmneplas(a);
   r_new(a,lmn)=lmnR(a);
   eta_new(a,lmn)=lmnETAP(a);
  end    
 end
end

if (graph==1 && step==SOL(1)) %Plot stress contours
 tri = delaunay(stress_pnt(:,1),stress_pnt(:,2));    
 figure
 hold on
 for i=1:size(tri,1)
  fill(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:)));
 end
 title('XFEM Stress')
 shading interp
 colorbar 
end
