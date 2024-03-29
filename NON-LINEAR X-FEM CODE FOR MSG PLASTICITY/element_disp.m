% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function U = element_disp(e,pos,enrich_node,u,element)

sctr=element(e,:);
nn=length(sctr);

idx=0;
stdU=zeros(2*nn,1);
for in=1:nn
 idx=idx+1;
 nodeI=sctr(in);
 stdU(2*idx-1)=u(2*nodeI-1);
 stdU(2*idx)=u(2*nodeI);
end

A=[];

if (any(enrich_node(sctr))==0)
 U=stdU;
else                             
 for in=1:nn
  nodeI=sctr(in);
  if (enrich_node(nodeI)==2)
   AA=[u(2*pos(nodeI)-1);u(2*pos(nodeI))];
   A=[A;AA];
  elseif(enrich_node(nodeI)==3)
   AA=[u(2*pos(nodeI)-1);u(2*pos(nodeI))];
   A=[A;AA];
  elseif(enrich_node(nodeI)==1)
   AA=[u(2*pos(nodeI)-1);u(2*pos(nodeI));u(2*(pos(nodeI)+1)-1);...
       u(2*(pos(nodeI)+1));u(2*(pos(nodeI)+2)-1);u(2*(pos(nodeI)+2));...
       u(2*(pos(nodeI)+3)-1);u(2*(pos(nodeI)+3));];
   A=[A;AA];
  end
 end
end

U = [stdU;A];
