% Function to compute the quadrature weights and locations
% Based on the original function by Jack Chessa

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [W,Q]=quadrature(IntOrder,intType,sdim)

if(nargin<3)
 if(strcmp(intType,'GAUSS')==1)
  sdim=1;
 else
  sdim = 2;
 end
end
if(nargin<2)
 intType='GAUSS';
end
if(strcmp(intType,'GAUSS')==1) 
    
 quadpoint=zeros(IntOrder^sdim ,sdim);
 quadweight=zeros(IntOrder^sdim,1);
  
 r1pt=zeros(IntOrder,1); 
 r1wt=zeros(IntOrder,1);

switch (IntOrder) 
 case 1
  r1pt(1)=0.000000000000000;
  r1wt(1)=2.000000000000000;

 case 2
  r1pt(1) = 0.577350269189626;
  r1pt(2) =-0.577350269189626;
  r1wt(1) = 1.000000000000000; 
  r1wt(2) = 1.000000000000000;         

 case 3
  r1pt(1) = 0.774596669241483;
  r1pt(2) =-0.774596669241483;
  r1pt(3) = 0.000000000000000;
  r1wt(1) = 0.555555555555556;
  r1wt(2) = 0.555555555555556; 
  r1wt(3) = 0.888888888888889;   

 case 4
  r1pt(1) = 0.861134311594053;
  r1pt(2) =-0.861134311594053;
  r1pt(3) = 0.339981043584856;
  r1pt(4) =-0.339981043584856;
  r1wt(1) = 0.347854845137454;
  r1wt(2) = 0.347854845137454; 
  r1wt(3) = 0.652145154862546;
  r1wt(4) = 0.652145154862546;  

 case 5
  r1pt(1) = 0.906179845938664;
  r1pt(2) =-0.906179845938664;
  r1pt(3) = 0.538469310105683;
  r1pt(4) =-0.538469310105683;
  r1pt(5) = 0.000000000000000;
  r1wt(1) = 0.236926885056189;
  r1wt(2) = 0.236926885056189;
  r1wt(3) = 0.478628670499366;
  r1wt(4) = 0.478628670499366;  
  r1wt(5) = 0.568888888888889;  

 case 6
  r1pt(1) = 0.932469514203152;
  r1pt(2) =-0.932469514203152;
  r1pt(3) = 0.661209386466265;
  r1pt(4) =-0.661209386466265;
  r1pt(5) = 0.238619186003152;
  r1pt(6) =-0.238619186003152;

  r1wt(1) = 0.171324492379170;
  r1wt(2) = 0.171324492379170;
  r1wt(3) = 0.360761573048139;
  r1wt(4) = 0.360761573048139;   
  r1wt(5) = 0.467913934572691; 
  r1wt(6) = 0.467913934572691;

 case 7
  r1pt(1) =  0.949107912342759;
  r1pt(2) = -0.949107912342759;
  r1pt(3) =  0.741531185599394;
  r1pt(4) = -0.741531185599394;
  r1pt(5) =  0.405845151377397;
  r1pt(6) = -0.405845151377397;
  r1pt(7) =  0.000000000000000;

  r1wt(1) = 0.129484966168870;
  r1wt(2) = 0.129484966168870;
  r1wt(3) = 0.279705391489277;
  r1wt(4) = 0.279705391489277;
  r1wt(5) = 0.381830050505119;
  r1wt(6) = 0.381830050505119;
  r1wt(7) = 0.417959183673469;

 case 8
  r1pt(1) =  0.960289856497536;
  r1pt(2) = -0.960289856497536;
  r1pt(3) =  0.796666477413627;
  r1pt(4) = -0.796666477413627;
  r1pt(5) =  0.525532409916329;
  r1pt(6) = -0.525532409916329;
  r1pt(7) =  0.183434642495650;
  r1pt(8) = -0.183434642495650;

  r1wt(1) = 0.101228536290376;
  r1wt(2) = 0.101228536290376;
  r1wt(3) = 0.222381034453374;
  r1wt(4) = 0.222381034453374;
  r1wt(5) = 0.313706645877887;
  r1wt(6) = 0.313706645877887;
  r1wt(7) = 0.362683783378362;
  r1wt(8) = 0.362683783378362;
   
 otherwise
  [W,Q] = disBlendingQ4quad(IntOrder,nodes); 
  return
 end  % end of IntOrder switch

 n=1;
     
 if(sdim==1) 
  for i=1:IntOrder
   quadpoint(n,:)=[r1pt(i)];           
   quadweight(n)=r1wt(i); 
   n=n+1;
  end
 elseif(sdim==2) 
  for i=1:IntOrder
   for j=1:IntOrder
    quadpoint(n,:)=[r1pt(i),r1pt(j)];           
    quadweight(n)=r1wt(i)*r1wt(j); 
    n=n+1;
   end
  end
 else % sdim == 3
  for i=1:IntOrder
   for j=1:IntOrder
    for k=1:IntOrder
     quadpoint(n,:)=[r1pt(i),r1pt(j),r1pt(k)];           
     quadweight(n)=r1wt(i)*r1wt(j)*r1wt(k); 
     n=n+1;
    end
   end
  end
      
 end
    
 Q=quadpoint;
 W=quadweight;
  
elseif (strcmp(intType,'TRIANGULAR')==1) 
    
if (sdim==3)  %%% TETRAHEDRA
      
 if (IntOrder~=1 && IntOrder~=2 && IntOrder~=3) 
  % check for valid quadrature order
  disp('Incorect quadrature order for triangular quadrature');
  IntOrder = 1;
 end
      
 if(IntOrder==1)
  quadpoint=[0.25 0.25 0.25];
  quadweight=1;
        
 elseif(IntOrder==2) 
  quadpoint=[0.58541020  0.13819660  0.13819660;
   0.13819660  0.58541020  0.13819660;
   0.13819660  0.13819660  0.58541020;
   0.13819660  0.13819660  0.13819660];
  quadweight = [1; 1; 1; 1]/4;
     
 elseif (IntOrder==3) 
quadpoint=[0.25 0.25 0.25;1/2 1/6 1/6;1/6 1/2 1/6;1/6 1/6 1/2;1/6 1/6 1/6];
  quadweight = [-4/5 9/20 9/20 9/20 9/20]';                      
 end
      
 Q=quadpoint;
 W=quadweight/6;
         
 else  %%% TRIANGLES
      
  if ( IntOrder > 7 ) % check for valid quadrature order
   disp('Quadrature order too high for triangular quadrature');
   IntOrder = 1;
  end
      
  if ( IntOrder == 1 )   % set quad points and quadweights
   quadpoint = [ 0.3333333333333, 0.3333333333333 ];
   quadweight = 1;
        
  elseif ( IntOrder == 2 ) 
   quadpoint = zeros( 3, 2 );
   quadweight = zeros( 3, 1 );
   
   quadpoint(1,:) = [ 0.1666666666667, 0.1666666666667 ];
   quadpoint(2,:) = [ 0.6666666666667, 0.1666666666667 ];
   quadpoint(3,:) = [ 0.1666666666667, 0.6666666666667 ]; 
   
   quadweight(1) = 0.3333333333333; 
   quadweight(2) = 0.3333333333333; 
   quadweight(3) = 0.3333333333333;   
        
  elseif ( IntOrder <= 5 ) 
   quadpoint = zeros( 7, 2 );
   quadweight = zeros( 7, 1 );
   
   quadpoint(1,:) = [ 0.1012865073235, 0.1012865073235 ];
   quadpoint(2,:) = [ 0.7974269853531, 0.1012865073235 ];
   quadpoint(3,:) = [ 0.1012865073235, 0.7974269853531 ]; 
   quadpoint(4,:) = [ 0.4701420641051, 0.0597158717898 ];
   quadpoint(5,:) = [ 0.4701420641051, 0.4701420641051 ];
   quadpoint(6,:) = [ 0.0597158717898, 0.4701420641051 ]; 
   quadpoint(7,:) = [ 0.3333333333333, 0.3333333333333 ];
   
   quadweight(1) = 0.1259391805448; 
   quadweight(2) = 0.1259391805448; 
   quadweight(3) = 0.1259391805448; 
   quadweight(4) = 0.1323941527885;
   quadweight(5) = 0.1323941527885;
   quadweight(6) = 0.1323941527885;
   quadweight(7) = 0.2250000000000;  
        
  else
   quadpoint = zeros( 13, 2 );
   quadweight = zeros( 13, 1 );
   

   quadpoint(1 ,:) = [ 0.0651301029022, 0.0651301029022 ];
   quadpoint(2 ,:) = [ 0.8697397941956, 0.0651301029022 ];
   quadpoint(3 ,:) = [ 0.0651301029022, 0.8697397941956 ];
   quadpoint(4 ,:) = [ 0.3128654960049, 0.0486903154253 ];
   quadpoint(5 ,:) = [ 0.6384441885698, 0.3128654960049 ];
   quadpoint(6 ,:) = [ 0.0486903154253, 0.6384441885698 ];
   quadpoint(7 ,:) = [ 0.6384441885698, 0.0486903154253 ];
   quadpoint(8 ,:) = [ 0.3128654960049, 0.6384441885698 ];
   quadpoint(9 ,:) = [ 0.0486903154253, 0.3128654960049 ];
   quadpoint(10,:) = [ 0.2603459660790, 0.2603459660790 ];
   quadpoint(11,:) = [ 0.4793080678419, 0.2603459660790 ];
   quadpoint(12,:) = [ 0.2603459660790, 0.4793080678419 ];
   quadpoint(13,:) = [ 0.3333333333333, 0.3333333333333 ];
   
   quadweight(1 ) = 0.0533472356088;
   quadweight(2 ) = 0.0533472356088; 
   quadweight(3 ) = 0.0533472356088;
   quadweight(4 ) = 0.0771137608903;
   quadweight(5 ) = 0.0771137608903;
   quadweight(6 ) = 0.0771137608903;
   quadweight(7 ) = 0.0771137608903;
   quadweight(8 ) = 0.0771137608903;
   quadweight(9 ) = 0.0771137608903;
   quadweight(10) = 0.1756152576332; 
   quadweight(11) = 0.1756152576332; 
   quadweight(12) = 0.1756152576332;
   quadweight(13) =-0.1495700444677; 
  end
      
  Q=quadpoint;
  W=quadweight/2;
 end
    
elseif ( strcmp(intType,'DUNAVANT') == 1 )
        
 if ( sdim == 3 )  %%% TETRAHEDRA
  disp ('not yet coded');
 else        
  [quadpoint,quadweight]=dunavant_rule(IntOrder );          
 end
   
 Q=quadpoint;
 W=quadweight/2;
end
