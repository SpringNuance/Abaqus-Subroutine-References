% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [W,Q] = disBlendQ4(IntOrder,subTriDiv,intType)

nodes=[-1 -1;1 -1;1 1;-1 1];

x=nodes(:,1);
y=nodes(:,2);

[geom]=polygeom(x,y);

node=[nodes;geom(2) geom(3)];

tri=delaunay(node(:,1),node(:,2)) ;
tri=tricheck(node,tri) ;

if(subTriDiv == 0)
 triangles=tri;
 triNodes=node;
else
 [triNodes,triangles]=subTriXFEM(node,tri,subTriDiv);
end

node = triNodes ;
tri = triangles ;

pt = 1;
for e=1:size(tri,1)
 [w,q]=quadrature(IntOrder,intType,2);
 coord=node(tri(e,:),:);
 a=det([coord,[1;1;1]])/2;
 if (a<0)
  coord=[coord(2,:);coord(1,:);coord(3,:)];
  a=det([coord,[1;1;1]])/2;
 end
 if(a~=0)
  for n=1:length(w)
   N=lagrange_basis('T3',q(n,:));
   Q(pt,:)=N'*coord;
   W(pt,1)=2*w(n)*a;
   pt=pt+1;
  end
 end
end