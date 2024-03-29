% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [W,Q] = disTipQ4(IntOrder,phi,nnode,tip,subTriDiv,elem,intType)

epsilon=0.00001;
corner=[1 2 3 4 1];
node=[-1 -1; 1 -1; 1 1; -1 1];

coord=zeros(1,2);
ksi=0;
eta=0;
iter=10;

inc=1;
while (inc<iter)
 [N,dNdxi]=lagrange_basis(elem,coord); 
 x=N'*nnode(:,1);
 y=N'*nnode(:,2);
 df1dr=dNdxi(:,1)'*nnode(:,1);
 df1ds=dNdxi(:,2)'*nnode(:,1);
 df2dr=dNdxi(:,1)'*nnode(:,2);
 df2ds=dNdxi(:,2)'*nnode(:,2);
 
 f1=x-tip(1);
 f2=y-tip(2);

 detF=df1dr*df2ds-df1ds*df2dr;

 invf(1,1)=1.0/detF*df2ds;
 invf(1,2)=-1.0/detF*df1ds;
 invf(2,1)=-1.0/detF*df2dr;
 invf(2,2)=1.0/detF*df1dr;

 ksi=ksi-invf(1,1)*f1-invf(1,2)*f2;
 eta=eta-invf(2,1)*f1-invf(2,2)*f2;

 coord(1)=ksi;
 coord(2)=eta;

 if((abs(ksi-coord(1))<epsilon)&&(abs(eta-coord(2))<epsilon))
  inc=iter+1;
  ntip=coord;
 else
  inc=inc+1;
 end
end

for i=1:4
 n1=corner(i);
 n2=corner(i+1);
 if phi(n1)*phi(n2) < 0
  r=phi(n1)/(phi(n1)-phi(n2));
  pnt=(1-r)*node(n1,:)+r*node(n2,:);
  node=[node;pnt];
 end
end

node=[node;ntip];

tri=delaunay(node(:,1),node(:,2)) ;
tri=tricheck(node,tri) ;

if(subTriDiv == 0)
 triangles=tri;
 triNodes=node;
else
 [triNodes,triangles]=subTriXFEM(node,tri,subTriDiv);
end

node=triNodes;
tri=triangles;

pt=1;
for e=1:size(tri,1)
 [w,q]=quadrature(IntOrder,intType,2);
 t=size(w,1);
 coord=node(tri(e,:),:);
 a=det([coord,[1;1;1]])/2;
 if(a<0)
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