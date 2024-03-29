function [B]=xfemBmat(pt,e,type_elem,enrich_node,xCrl,GVertex,cont,node,...
    element,MAT,tip_elem,phiN,elem,blend)

sctr=element(e,:);
nn=length(sctr);
[N,dNdxi]=lagrange_basis(elem,pt);
J0=node(sctr,:)'*dNdxi;             
dNdx=dNdxi*inv(J0);    
Gpt=N'*node(sctr,:);

if blend==1
 rampFn=phiN(sctr,:);   
 Rgp=N'*rampFn ;
 dRgp=dNdx'*rampFn ;
else
 Rgp=1;
 dRgp=zeros(2,1);
end

if cont==1
 Bfem=zeros(3,2*nn);
 Bfem(1,1:2:2*nn)=dNdx(:,1)';
 Bfem(2,2:2:2*nn)=dNdx(:,2)';
 Bfem(3,1:2:2*nn)=dNdx(:,2)';
 Bfem(3,2:2:2*nn)=dNdx(:,1)';
else
 Bfem=[];
end

if(any(enrich_node(sctr))==0)
 B=Bfem ;
else
 Bxfem=[];
 for in=1:nn
  if(enrich_node(sctr(in))==2)
   if(type_elem(e,cont)==2) || (type_elem(e,cont)==3)
    e_ref=e;
    xCre=[xCrl(e_ref,1) xCrl(e_ref,2); xCrl(e_ref,3) xCrl(e_ref,4)];
    if(type_elem(e,cont)==3)
     distV=signed_distance(xCre,GVertex(e_ref,:),0);
     HV=sign(distV);
     dist=signed_distance(xCre,Gpt,0);
     Hgp=sign(dist);
     if HV*Hgp<=0
      Hgp=Hgp;
     else
      vv=[xCre(1,:);xCre(2,:);GVertex(e_ref,:)];
      flag=inhull(Gpt,vv);
      if flag==1
       Hgp=-Hgp;
      else    
       Hgp=Hgp;
      end
     end
     dist=signed_distance(xCre,node(sctr(in),:),0);
     Hi=sign(dist);               
    else 
     dist=signed_distance(xCre,Gpt,16);
     Hgp=sign(dist);
     dist=signed_distance(xCre,node(sctr(in),:),0);
     Hi=sign(dist);
    end
    at1=dNdx(in,1)*Rgp*(Hgp-Hi)+N(in)*dRgp(1)*(Hgp-Hi);
    at2=dNdx(in,2)*Rgp*(Hgp-Hi)+N(in)*dRgp(2)*(Hgp-Hi);
    BI_enr=[at1 0;0 at2;at2 at1];
   else
    BI_enr=[0 0; 0 0; 0 0];
   end
   Bxfem = [Bxfem BI_enr];
   clear BI_enr;          
  elseif(enrich_node(sctr(in))==1) 
   e_ref=tip_elem(1);                
   xCre=[xCrl(e_ref,1) xCrl(e_ref,2); xCrl(e_ref,3) xCrl(e_ref,4)];
   seg=xCre(2,:)-xCre(1,:);
   alpha=atan2(seg(2),seg(1));
   xTip=[xCre(2,1) xCre(2,2)];
   QT=[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
   xp=QT*(Gpt-xTip)';          
   r=sqrt(xp(1)*xp(1)+xp(2)*xp(2));
   theta=atan2(xp(2),xp(1));
   if(theta>pi || theta<-pi)
    disp(['something wrong with angle ',num2str(thet)]);
   end
   [Br,dBdx,dBdy]=branch_gp(r,theta,alpha,MAT);
   xp=QT*(node(sctr(in),:)-xTip)';
   r=sqrt(xp(1)*xp(1)+xp(2)*xp(2));
   theta=atan2(xp(2),xp(1));   
   if(theta>pi || theta<-pi)
    disp(['something wrong with angle ',num2str(thet)]);
   end
   [BrI]=branch_node(r,theta,MAT);
   aa=dNdx(in,1)*(Br(1)-BrI(1))*Rgp+N(in)*dBdx(1)*Rgp+...
       N(in)*(Br(1)-BrI(1))*dRgp(1);
   bb=dNdx(in,2)*(Br(1)-BrI(1))*Rgp+N(in)*dBdy(1)*Rgp+...
       N(in)*(Br(1)-BrI(1))*dRgp(2);  
   B1_enr=[aa 0; 0 bb; bb aa];
   aa=dNdx(in,1)*(Br(2)-BrI(2))*Rgp+N(in)*dBdx(2)*Rgp+...
       N(in)*(Br(2)-BrI(2))*dRgp(1);
   bb=dNdx(in,2)*(Br(2)-BrI(2))*Rgp + N(in)*dBdy(2)*Rgp+...
       N(in)*(Br(2)-BrI(2))*dRgp(2);
   B2_enr=[aa 0; 0 bb; bb aa];
   aa=dNdx(in,1)*(Br(3)-BrI(3))*Rgp+N(in)*dBdx(3)*Rgp+...
       N(in)*(Br(3)-BrI(3))*dRgp(1);
   bb=dNdx(in,2)*(Br(3)-BrI(3))*Rgp+N(in)*dBdy(3)*Rgp+...
       N(in)*(Br(3)-BrI(3))*dRgp(2);
   B3_enr=[aa 0; 0 bb; bb aa];
   aa=dNdx(in,1)*(Br(4)-BrI(4))*Rgp+N(in)*dBdx(4)*Rgp+...
       N(in)*(Br(4)-BrI(4))*dRgp(1);
   bb=dNdx(in,2)*(Br(4)-BrI(4))*Rgp+N(in)*dBdy(4)*Rgp+...
       N(in)*(Br(4)-BrI(4))*dRgp(2);
   B4_enr=[aa 0; 0 bb; bb aa];         
   BI_enr=[B1_enr B2_enr B3_enr B4_enr];
   clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
   Bxfem = [Bxfem BI_enr];
   clear BI_enr ;
  end
 end 
 B=[Bfem Bxfem];
 clear Bfem; clear Bxfem;
end