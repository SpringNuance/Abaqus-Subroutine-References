% Function to plot the finite element mesh

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function plotmesh(coords,connect,color,elem)

hold on

if(size(coords,2)<3)
 for i=size(coords,2)+1:3
  coords(:,i)=zeros(size(coords,1),1);
 end
end

for i=1:size(connect,1)
 if(strcmp(elem,'Q8'))
  order=[1,5,2,6,3,7,4,8,1];
 elseif(strcmp(elem,'Q4'))
  order=[1,2,3,4,1]; 
 end
 xpt=zeros(size(order,2));  
 ypt=zeros(size(order,2));
 zpt=zeros(size(order,2));
 for j=1:size(order,2)
  xpt(j)=coords(connect(i,order(j)),1);
  ypt(j)=coords(connect(i,order(j)),2);      
  zpt(j)=coords(connect(i,order(j)),3);
 end
 plot3(xpt,ypt,zpt,color)
end

axis equal