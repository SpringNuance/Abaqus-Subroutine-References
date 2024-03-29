% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function sctrB=assembly(lmn,enrich_node,pos,element)

sctr=element(lmn,:);
nn=length(sctr);

for k=1:nn
 sctrBfem(2*k-1)=2*sctr(k)-1;
 sctrBfem(2*k)=2*sctr(k);
end

if(any(enrich_node(sctr))==0) 
 sctrB=sctrBfem;
else
 tn=size(find(enrich_node(sctr)==1),1);
 sn=size(find(enrich_node(sctr)==2),1);
 ssn=size(find(enrich_node(sctr)==3),1);
 sctrBxfem=zeros(1,2*(sn*1+tn*4+ssn*1));
 cnt=0;
 for k=1:nn
  if(enrich_node(sctr(k))==2) 
   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*pos(sctr(k))-1;
   sctrBxfem(2*cnt)=2*pos(sctr(k));
  elseif(enrich_node(sctr(k))==3)  
   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*pos(sctr(k))-1;
   sctrBxfem(2*cnt)=2*pos(sctr(k));
  elseif(enrich_node(sctr(k))==1) 
   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*pos(sctr(k))-1;
   sctrBxfem(2*cnt)=2*pos(sctr(k));

   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*(pos(sctr(k))+1)-1;
   sctrBxfem(2*cnt)=2*(pos(sctr(k))+1);

   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*(pos(sctr(k))+2)-1;
   sctrBxfem(2*cnt)=2*(pos(sctr(k))+2);

   cnt=cnt+1;
   sctrBxfem(2*cnt-1)=2*(pos(sctr(k))+3)-1;
   sctrBxfem(2*cnt)=2*(pos(sctr(k))+3);
  end
 end
 sctrB=[sctrBfem sctrBxfem];
end
