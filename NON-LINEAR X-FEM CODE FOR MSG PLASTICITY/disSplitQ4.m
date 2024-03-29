% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [W,Q]=disSplitQ4(IntOrder,phi,subTriDiv,intType)

corner=[1 2 3 4 1];
node=[-1 -1;1 -1;1 1;-1 1];

numEdge=size(node,1);
cutEdge=[];
for i = 1:numEdge
 n1=corner(i);
 n2=corner(i+1);
 if(phi(n1)*phi(n2)<0)
  r=phi(n1)/(phi(n1)-phi(n2));
  pnt=(1-r)*node(n1,:)+r*node(n2,:);
  node=[node;pnt];
  cutEdge=[cutEdge i];
 end
end

nEdge=length(cutEdge);
if(cutEdge(2)==cutEdge(1)+1 || cutEdge(2)==cutEdge(1)+3)
 if(ismember(cutEdge(1),[1 2]) && ismember(cutEdge(2),[1,2]))
  tempNode=[node(1,:);node(5,:);node(6,:);node(3,:);node(4,:)];
  [geom]=polygeom(tempNode(:,1),tempNode(:,2));
  node=[node;geom(2) geom(3)];
  tri=[1 7 5;1 7 4;4 7 3;3 7 6;6 7 5;5 6 2];
  tri=tricheck(node,tri);
 elseif(ismember(cutEdge(1),[2 3]) && ismember(cutEdge(2),[2 3]))
  tempNode=[node(1,:);node(2,:);node(5,:);node(6,:);node(4,:)];
  [geom]=polygeom(tempNode(:,1),tempNode(:,2));
  node=[node;geom(2) geom(3)];
  tri=[1 7 2;1 7 4;4 7 6;6 7 5;5 7 2;3 6 5];
  tri=tricheck(node,tri);
 elseif(ismember(cutEdge(1),[3 4]) && ismember(cutEdge(2),[3 4]))
  tempNode=[node(1,:);node(2,:);node(3,:);node(5,:);node(6,:)];
  [geom]=polygeom(tempNode(:,1),tempNode(:,2));
  node=[node;geom(2) geom(3)];
  tri=[1 2 7;2 7 3;7 3 5;5 6 7;6 7 1;4 5 6];
  tri=tricheck(node,tri);
 elseif(ismember(cutEdge(1),[1 4]) && ismember(cutEdge(2),[1 4]))
  tempNode=[node(5,:);node(2,:);node(3,:);node(4,:);node(6,:)];
  [geom]=polygeom(tempNode(:,1),tempNode(:,2));
  node=[node;geom(2) geom(3)];
  tri=[1 5 6;5 7 2;6 7 5;6 7 4;4 7 3;7 3 2];
  tri=tricheck(node,tri);
 end
else
 if(cutEdge(1)==2 && cutEdge(2)==4)
  tempNode1=[node(1,:);node(2,:);node(5,:);node(6,:)];
  tempNode2=[node(6,:);node(5,:);node(3,:);node(4,:)];

  [geom]=polygeom(tempNode1(:,1),tempNode1(:,2));
  CenPoint1=[geom(2) geom(3)];

  [geom]=polygeom(tempNode2(:,1),tempNode2(:,2));
  CenPoint2=[geom(2) geom(3)];
  node=[node; CenPoint1; CenPoint2];
  tri=[1 7 2;2 5 7;7 5 6;6 7 1;6 8 5;5 8 3;3 8 4;4 8 6];
  tri=tricheck(node,tri);
 elseif(cutEdge(1)==1 && cutEdge(2)==3)
  tempNode1=[node(1,:);node(5,:);node(6,:);node(4,:)];
  tempNode2=[node(5,:);node(2,:);node(3,:);node(6,:)];
  [geom]=polygeom(tempNode1(:,1),tempNode1(:,2));
  CenPoint1=[geom(2) geom(3)];

  [geom]=polygeom(tempNode2(:,1),tempNode2(:,2));
  CenPoint2=[geom(2) geom(3)];
  node=[node; CenPoint1; CenPoint2];

  tri=[1 7 4;1 7 5;7 5 6;4 7 6;6 8 5;5 8 2;8 2 3;3 8 6];
  tri=tricheck(node,tri);
 end
end

if(subTriDiv==0)
 triangles=tri;
 triNodes=node;
else
 [triNodes,triangles]=subTriXFEM(node,tri,subTriDiv);
end

node=triNodes;
tri=triangles;

pt = 1;
for i=1:size(tri,1)
 [w,q]=quadrature(IntOrder,intType,2);
 coord=node(tri(i,:),:);
 a=det([coord,[1;1;1]])/2;
 if(a<0)
  coord=[coord(2,:);coord(1,:);coord(3,:)];
  a=det([coord,[1;1;1]])/2;
 end
 if(a~=0)
  for n=1:length(w)
   N=lagrange_basis('T3',q(n,:));
   Q(pt,:) = N'*coord;
   W(pt,1) = 2*w(n)*a;
   pt = pt+1;
  end
 end
end
