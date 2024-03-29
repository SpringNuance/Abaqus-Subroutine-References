% Function to store integration points and weights for each element

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [W,Q]=gauss_rule(e,enrich_node,elem_crk,xTip,xVertex,tip_elem,...
    split_elem,vertex_elem,node,element,elem)

sctr=element(e,:);
nnode=node(sctr,:);

tip_enr=find(enrich_node(sctr,:)==1);
tipnodes=[];
tipnodes=[tipnodes; find(enrich_node(:,1)==1)];

EPS=1e-12;
x0=elem_crk(e,1); 
y0=elem_crk(e,2);
x1=elem_crk(e,3); 
y1=elem_crk(e,4);

for i=1:size(sctr,2)
 x=node(sctr(i),1);
 y=node(sctr(i),2);
 phiT=(y0-y1)*x+(x1-x0)*y+(x0*y1-y0*x1);
 if abs(phiT)<EPS
  phi(i,1)=0;
 else
  phi(i,1)=phiT;
 end
end

NormalOrder=2;  
TipOrder=6;     
SplitOrder=3;
VertexOrder=3;

if(ismember(e,tip_elem))
 IntOrder=TipOrder;
elseif(ismember(e,split_elem) && any(ismember(sctr,tipnodes)))
 IntOrder=TipOrder;
elseif(ismember(e,vertex_elem))
 IntOrder=VertexOrder;
elseif(size(tip_enr,1)>0)
 IntOrder=TipOrder;
elseif( ismember(e,split_elem) )
 for c=1:length(split_elem)
  M=find(e==split_elem(c));
  if(M==1)
   IntOrder=SplitOrder;
  end
 end
else
 IntOrder=NormalOrder;
end

if(IntOrder==TipOrder && TipOrder<=7)
 intType='DUNAVANT';
 subTriDiv=0;
elseif(IntOrder==TipOrder && TipOrder>7)
 intType='DUNAVANT';
 subTriDiv=0;
end

if(ismember(e,tip_elem))
 [W,Q]=disTipQ4(IntOrder,phi,nnode,xTip(e,:),subTriDiv,elem,intType);
elseif(ismember(e,split_elem) && any(ismember(sctr,tipnodes)))
 [W,Q]=disSplitQ4(IntOrder,phi,subTriDiv,intType);
elseif(ismember(e,vertex_elem))
  subTriDiv=0;
  intType='TRIANGULAR';
  [W,Q]=disTipQ4(IntOrder,phi,nnode,xVertex(e,:),subTriDiv,elem,intType);
elseif(size(tip_enr,1)>0)
  [W,Q]=disBlendQ4(IntOrder,subTriDiv,intType);
elseif(ismember(e,split_elem))
 for c=1:length(split_elem)
  M=find(e==split_elem(c));
  if(M==1)
   subTriDiv=0;
   intType='TRIANGULAR';
   [W,Q]=disSplitQ4(IntOrder,phi,subTriDiv,intType);
  end
 end
else
 intType='GAUSS';
 [W,Q]=quadrature(IntOrder,intType,2);
end