% Function to compute the intersection of two line segments in 2D.
% Based on the original function by John Burkardt

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function flag=segments_int_2d(p1,p2,q1,q2)

dim_num=2;
tol=1e-7; %IMPORTANT (mesh-dependent): decrease with the tip element size
r=[];

ival=0;
r(1:dim_num)=0;

%  Check whether either line is a point.
if(p1(1:dim_num)==p2(1:dim_num))
 point_1=1;
else
 point_1=0;
end
if(q1(1:dim_num)==q2(1:dim_num))
 point_2=1;
else
 point_2=0;
end

%  Convert the lines to ABC format.
if(~point_1)
 if (p1(1:dim_num)==p2(1:dim_num))
  fprintf(1,'\n');
  fprintf(1,'LINE_EXP2IMP_2D - Fatal error!\n');
  fprintf(1,'P1 = P2\n');
  fprintf(1,'P1 = %f %f\n',p1(1:dim_num));
  fprintf(1,'P2 = %f %f\n',p2(1:dim_num));
  error ('LINE_EXP2IMP_2D - Fatal error!');
 end

 a1=p2(2)-p1(2);
 b1=p1(1)-p2(1);
 c1=p2(1)*p1(2)-p1(1)*p2(2);

 norm=a1*a1+b1*b1+c1*c1;
 if(0<norm)
  a1=a1/norm;
  b1=b1/norm;
  c1=c1/norm;
 end

 if(a1<0)
  a1=-a1;
  b1=-b1;
  c1=-c1;
 end
end

if(~point_2)  
 if (q1(1:dim_num)==q2(1:dim_num))
  fprintf(1,'\n');
  fprintf(1,'LINE_EXP2IMP_2D - Fatal error!\n');
  fprintf(1,'Q1 = Q2\n');
  fprintf(1,'Q1 = %f %f\n',q1(1:dim_num));
  fprintf(1,'Q2 = %f %f\n',q2(1:dim_num));
  error ('LINE_EXP2IMP_2D - Fatal error!');
 end

 a2=q2(2)-q1(2);
 b2=q1(1)-q2(1);
 c2=q2(1)*q1(2)-q1(1)*q2(2);

 norm=a2*a2+b2*b2+c2*c2;
 if(0<norm)
  a2=a2/norm;
  b2=b2/norm;
  c2=c2/norm;
 end

 if(a2<0)
  a2=-a2;
  b2=-b2;
  c2=-c2;
 end 
end

%  Search for intersection of the lines.
if(point_1 && point_2)
 if(p1(1:dim_num)==q1(1:dim_num))
  ival=1;
  r(1:dim_num)=p1(1:dim_num);
 end
elseif(point_1)
 if(a2*p1(1)+b2*p1(2)==c2)
  ival=1;
  r(1:dim_num)=p1(1:dim_num);
 end
elseif(point_2)
 if(a1*q1(1)+b1*q1(2)==c1)
  ival=1;
  r(1:dim_num)=q1(1:dim_num);
 end
else
 if(a1==0 && b1==0)
  ival=-1;
  r=[];
 elseif(a2==0 && b2==0)
  ival=-2;
  r=[];
 else
  a(1,1)=a1;
  a(1,2)=b1;
  a(1,3)=-c1;

  a(2,1)=a2;
  a(2,2)=b2;
  a(2,3)=-c2;
  
  n=2;
  nrhs=1;
  info=0;

  for j=1:n
%  Choose a pivot row IPIVOT.
   ipivot = j;
   apivot = a(j,j);

   for i=j+1:n
    if (abs(apivot)<abs(a(i,j)))
     apivot=a(i,j);
     ipivot=i;
    end
   end

   if(apivot==0.0)
    info = j;
    break
   else
%  Interchange.
    temp=a(ipivot,1:n+nrhs);
    a(ipivot,1:n+nrhs)=a(j,1:n+nrhs);
    a(j,1:n+nrhs)=temp;
%  A(J,J) becomes 1.
    a(j,j)=1.0;
    a(j,j+1:n+nrhs)=a(j,j+1:n+nrhs)/apivot;
%  A(I,J) becomes 0.
    for i=1:n
     if(i~=j)
      factor=a(i,j);
      a(i,j)=0.0;
      a(i,j+1:n+nrhs)=a(i,j+1:n+nrhs)-factor*a(j,j+1:n+nrhs);
     end
    end
   end       
  end
  
% If the inverse exists, then the lines intersect at the solution point.
  if(info==0)
   ival=1;
   r(1:dim_num)=a(1:dim_num,3);
% If the inverse does not exist, then the lines are parallel or coincident.  
% Check for parallelism by seeing if the C entries are in the same ratio 
% as the A or B entries.   
  else
   ival=0;
   r=[];
   if(a1==0)
    if(b2*c1==c2*b1)
     ival=2;
     r(1:dim_num)=[0,-c1/b1];
    end
   else
    if(a2*c1==c2*a1)
     ival=2;
     if(abs(a1)<abs(b1))
      r(1:dim_num)=[0,-c1/b1];
     else
      r(1:dim_num)=[-c1/a1,0];
     end
    end
   end
  end
 end 
end

if (ival==0)
 flag = 0 ;
 return
end

% Is the point on the first segment?
normsq=sum((p2(1:dim_num)-p1(1:dim_num)).^2);
if(normsq==0.0)
 if(r(1:dim_num)==p1(1:dim_num))
  u=0.5;
  v=0.0;
 else
  u=0.5;
  v=Inf;
 end
else
 u=((r(1:dim_num)-p1(1:dim_num))*(p2(1:dim_num)-p1(1:dim_num))')/normsq;
 v=sqrt(((u-1.0)*p1(1)-u*p2(1)+r(1)).^2+((u-1.0)*p1(2)-u*p2(2)+r(2)).^2)...
   /sqrt(normsq);
end

if(u <0.0-tol||1.0+tol<u||tol<v)
 flag=0;
 return
end

% Is the point on the second segment?
normsq=sum((q2(1:dim_num)-q1(1:dim_num)).^2);
if(normsq==0.0)
 if(r(1:dim_num)==q1(1:dim_num))
  u=0.5;
  v=0.0;
 else
  u=0.5;
  v=Inf;
 end
else
 u=((r(1:dim_num)-q1(1:dim_num))*(q2(1:dim_num)-q1(1:dim_num))')/normsq;
 v=sqrt(((u-1.0)*q1(1)-u*q2(1)+r(1)).^2+((u-1.0)*q1(2)-u*q2(2)+r(2)).^2)...
   /sqrt(normsq);  
end  

if(u<0.0-tol||1.0+tol<u||tol<v)
 flag = 0;
 return
end

flag = [1 r];