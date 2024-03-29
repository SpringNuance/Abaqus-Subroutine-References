%Function to compute the residual and the global stiffness matrix

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [GK,resid,stress_pnt,STRAINP,strainp_val,F1]=globalKR(dt,ndof,...
 coords,nelem,connect,materialprops,stress,eplas,w,MAT,ee,ep,eta,nne,tu,...
 enrich_node,elem_crk,type_elem,xTip,xVertex,split_elem,tip_elem,pos,...
 vertex_elem,step,nit,phiN,elem,blend,stress_pnt,STRAINP,strainp_val,F1)

node=coords';
element=connect';
[NCorN, ~]=size(node);

% Assemble the global stiffness matrix
resid=zeros(tu,1);
IJX=zeros(ndof*nne*nelem*ndof*nne,3);
count=1;
lmncoord=zeros(2,nne);
lmndof=zeros(ndof,nne);
   
% Empty strainp1_val, make strainp1_valTemp=strainp1
STRAINPT=zeros(NCorN,4);
Icount=1;
   
if MAT==2
 if step==1 && nit==2
  test=zeros(size(stress_pnt,1),1);
  F1=scatteredInterpolant(stress_pnt(:,1),stress_pnt(:,2),test);
 end     
 if step>1 || nit>2
  F1.Values=strainp_val(:,1);
  STRAINPT(:,1) = F1(node);
  F1.Values=strainp_val(:,2); 
  STRAINPT(:,2) = F1(node);
  F1.Values=strainp_val(:,3);
  STRAINPT(:,3) = F1(node);
  F1.Values=strainp_val(:,4);
  STRAINPT(:,4) = F1(node);
  if step~=1 && nit~=1
   strainp_val=zeros(size(stress_pnt,1),4);
  end    
 end
end

for lmn=1:nelem
       
 sctr = element(lmn,:);  
 coordN=node(sctr',:);
 [Nlines, ~] = size(coordN);
 if Nlines>4
  coordN(5,:)=[];coordN(5,:)=[];coordN(5,:)=[];coordN(5,:)=[];
 end   
   
 if MAT==2
  if step>1 || nit>2
   sctrT=sctr(1,1:4);
   STRAINP(lmn,:,1) = STRAINPT(sctrT',1);
   STRAINP(lmn,:,2) = STRAINPT(sctrT',2);
   STRAINP(lmn,:,3) = STRAINPT(sctrT',3);  
   STRAINP(lmn,:,4) = STRAINPT(sctrT',4);
   end
  end   
   
 [W,Q]=gauss_rule(lmn,enrich_node,elem_crk,xTip,xVertex,tip_elem,...
     split_elem,vertex_elem,node,element,elem);   
 
 for igp=1:size(W,1)
  gpnt=Q(igp,:) ;
  [N,~]=lagrange_basis(elem,gpnt) ;
   pt=N'*node(sctr,:);
        
  if step==1 && nit==1
   stress_pnt = [stress_pnt; pt] ; 
  end
 end   
      
 sctrB=[];
 sctrB=[sctrB assembly(lmn,enrich_node(:,1),pos(:,1),element)];   
        
 U=[];
 U=[U; element_disp(lmn,pos(:,1),enrich_node(:,1),w,element)];  
    
 for a=1:nne
  for i=1:2
   lmncoord(i,a)=coords(i,connect(a,lmn));
  end
  for i=1:ndof
   lmndof(i,a)= w(ndof*(connect(a,lmn)-1)+i);
  end
 end
 nintp=size(W,1);
 lmnstress=zeros(4,nintp);
 lmneplas=zeros(nintp,1);
 lmnEPL=zeros(4,nintp);
 lmnEE=zeros(4,nintp);
 lmnETAP=zeros(nintp,1);
 if MAT > 0
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
  
 Kel=zeros(size(sctrB,2),size(sctrB,2)); 
 Rel=zeros(size(sctrB,2),1);

 for kk=1:size(W,1)
  B=[];
  Gpt=Q(kk,:);
  [~,dNdxi]=lagrange_basis(elem,Gpt);
  JO=node(sctr,:)'*dNdxi;
  B=[B xfemBmat(Gpt,lmn,type_elem,enrich_node(:,1),elem_crk,xVertex,1,...
    node,element,MAT,tip_elem,phiN,elem,blend)] ;
  eps_sub=B*U ;
     
  [C,stresst,strainp_val,Icount]=ElementKR(dt,ndof,materialprops,...
      lmnstress,lmneplas,MAT,lmnETAP,lmnEE,lmnEPL,lmn,kk,eps_sub,coordN,...
      Gpt,STRAINP,strainp_val,Icount);
          
  Kel=Kel+B'*C*B*W(kk)*det(JO);
  Rel=Rel+B'*stresst'*W(kk)*det(JO);
 end %Gauss points
 B = repmat(sctrB,size(sctrB,2),1);
 C=reshape(B,[1],[])';
 D=reshape(B',[],[1]);
 Kel=reshape(Kel,[],[1]);
 IJX(count:(count+size(sctrB,2)^2-1),:)=horzcat(C,D,Kel);
 count=count+size(sctrB,2)^2;
 resid(sctrB) = resid(sctrB) + Rel;
end %Elements
IJX = IJX(all(IJX,2),:);
I=IJX(:,1);J=IJX(:,2);X=IJX(:,3);
GK=sparse(I,J,X,tu,tu);