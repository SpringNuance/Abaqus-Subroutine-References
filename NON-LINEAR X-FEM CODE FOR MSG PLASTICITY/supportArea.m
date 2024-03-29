% Function to test the tolerance linked to the support area

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [Aw,Awp]=supportArea(pt,e,TypeEl,elem_crk,xVertex,node,elem)

[sctrn,~]=find(elem == pt);      % find support elements
coor_pt=node(pt,:);                      
xCre=[elem_crk(e,1) elem_crk(e,2); elem_crk(e,3) elem_crk(e,4)];
distPT=signed_distance(xCre,coor_pt,1) ;  % below or above?  (sign)
HPT=sign(distPT)  ;                     % below or above?  (-1 or 1)

Awp=0; % area +
Aw=0; % total area

for i=1:size(sctrn,1)                                    
 sctr=elem(i,:);                                     
 Aw=Aw+polyarea(node(sctr,1),node(sctr,2));            
end               

if (HPT == 0) || (abs(distPT) < 1e-4)  
 HPT=1;
end

for i = 1: size(sctrn,1)           
 sctr=elem(sctrn(i),:);
 nn=length(sctr);
 pts=[];
     
 if (TypeEl(sctrn(i),1) == 0)
  for in = 1:nn                    
   nd=node(sctr(in),:);
   dist=signed_distance(xCre,nd,1);
   Hi=sign(dist);
   if Hi==HPT                
    pts = [pts; nd];
   end
  end
 end
 
 if (TypeEl(sctrn(i),1) == 2) || (TypeEl(sctrn(i),1) == 3) % cut
  xCr_local=[elem_crk(sctrn(i),1) elem_crk(sctrn(i),2);...
                      elem_crk(sctrn(i),3) elem_crk(sctrn(i),4)];   
  for in = 1:nn  
   nd = node(sctr(in),:);
   dist = signed_distance(xCr_local,nd,1);
   Hi = sign(dist);
   if Hi == HPT 
    pts = [pts; nd];
   end
  end
  pts = [pts; xCr_local];  
 end
 
 if (size(pts,1) > 2 )              
  [~,surf] = convhull(pts(:,1),pts(:,2));  
  Awp = Awp + surf;
 else    
  Awp = Awp;
 end
    
 if (TypeEl(sctrn(i),1) == 3) 
  dist = signed_distance(xCr_local,xVertex(sctrn(i),:),1);
  Hi = sign(dist);
  points = [xCr_local ; xVertex(sctrn(i),:)];
  [~,surf] = convhull(points(:,1),points(:,2));
  Awp=Awp-Hi*HPT*surf;   
 end       
end

