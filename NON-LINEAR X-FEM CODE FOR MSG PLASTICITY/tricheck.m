% Function to check if a triangle has a negative Jacobian

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function conn=tricheck(node,conn,verbose)

if(nargin==2)
 verbose=0;
end
if(size(node,2)==3)
 node=node(:,1:2);
end

count=0;

for e=1:size(conn,1)
 sctr=conn(e,:);
 [~,dNdxi]=lagrange_basis('T3',[1/3 1/3]);
 detJ=det(node(sctr,:)'*dNdxi);
  
 if (detJ<0)
  conn(e,:)=fliplr(sctr);
  count=count+1;
 elseif ( detJ == 0 )
  disp(['ZERO JACOBIAN IN ELEMENT ',num2str(e),' CANNOT FIX'])
 end
end

if(verbose)
 disp(['TRICHECK FOUND ',num2str(count),' NEGATIVE JACOBIANS, ALL FIXED'])
end