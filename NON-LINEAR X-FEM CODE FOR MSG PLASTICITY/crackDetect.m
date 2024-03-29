% Function to identify the enriched domain in the vicinity of the crack tip

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [typeEl,ElCrk,tipEl,splitEl,vertexEl,xTip,xVertex,enrNode]=...
    crackDetect(crack,XY,LE)

Tip1=crack(1,:);
Tip2=crack(size(crack,1),:);
r1=[];
r2=[];
for i=1:size(XY,1)
 sctr=XY(i,:);
 rho1=sqrt((sctr(1)-Tip1(1))^2+(sctr(2)-Tip1(2))^2);
 r1=[r1,rho1];
 rho2=sqrt((sctr(1)-Tip2(1))^2+(sctr(2)-Tip2(2))^2);
 r2=[r2,rho2]; 
end
test1=r1-sqrt((Tip2(1)-Tip1(1))^2+(Tip2(2)-Tip1(2))^2);
test2=r2-sqrt((Tip1(1)-Tip2(1))^2+(Tip1(2)-Tip2(2))^2);
test1=test1(LE)';
test2=test2(LE)';
test1=min(test1);
test2=min(test2);
enrdomain1=find(test1<=0)';
enrdomain2=find(test2<=0)';
enrdomain=[enrdomain1;enrdomain2];
enrdomain=unique(enrdomain);

typeEl=zeros(size(LE,1),1);
ElCrk=zeros(size(LE,1),4);
xCrelement=zeros(size(LE,1),2);
xTip=zeros(size(LE,1),2);
xVertex=zeros(size(LE,1),2);
enrNode=zeros(size(XY,1),1);
tipEl=[];
splitEl=[];
vertexEl=[];

%Identify enriched elements(tip, vertex, split)
for i=1:size(enrdomain,1)
 e=enrdomain(i);
 sctr=LE(e,:);
 vv=XY(sctr,:);
 crkInt=[];  
 intes=0;
 flag1=0;
 flag2=0;
 for j=1:size(crack,1)-1       
  q1=crack(j,:); 
  q2=crack(j+1,:);
  sctrl=[sctr sctr(1,1)];      
  for k=1:size(sctr,2)
   nnode1=sctrl(k);
   nnode2=sctrl(k+1);
   p1 = XY(nnode1,:);
   p2 = XY(nnode2,:);
   intersect=segments_int_2d(p1,p2,q1,q2) ;
   intes = intes + intersect(1);
   if intersect(1) > 0
    crkInt = [crkInt intersect(2) intersect(3)];
    flag1 = inhull(crack(j,:),vv,[],-1e-8);
    flag2 = inhull(crack(j+1,:),vv,[],-1e-8);
    xCrelement(e,:)=crack(j,:)*flag1+crack(j+1,:)*flag2;
   end 
  end % k (edges of the elements)
 end % j (crack elements)
         
 if((intes==2) && (flag1==0) && (flag2==0))     % SPLIT
  typeEl(enrdomain(i))=2; 
  splitEl=[splitEl;e];
  ElCrk(e,:)=crkInt;
 end
 if (intes == 1)                                     % TIP
  typeEl(e)=1;  
  tipEl=[tipEl;e];
  xTip(e,:)=xCrelement(e,:);
  ElCrk(e,:)=[crkInt xTip(e,1) xTip(e,2)];      
 end
end % i (elements)   

% select the enriched nodes  
for i=1:size(enrdomain,1)
 sctr = LE(enrdomain(i),:);
 if typeEl(enrdomain(i))==1        % tip
  enrNode(sctr)=1;
 elseif typeEl(enrdomain(i))==2   % split
  for in=1:length(sctr)               % loop on the nodes of the element
   if enrNode(sctr(in))==0 % already enriched
    [Aw,Awp]=supportArea(sctr(in),enrdomain(i),typeEl,ElCrk,xVertex,XY,LE);
    if (abs(Awp / Aw) > 1e-4) && (abs((Aw-Awp) / Aw) > 1e-4)   
     enrNode(sctr(in))=2;
    end
   end
  end
 end  %if
end  % i (nodes)

% Emilio Martínez-Pañeda (mail@empaneda.com)
% www.empaneda.com/codes

