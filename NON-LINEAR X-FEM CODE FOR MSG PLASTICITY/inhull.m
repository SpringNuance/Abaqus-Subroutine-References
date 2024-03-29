% Function to test if a set of points is inside of a convex hull
% Based on the original function by John D'Errico

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function in = inhull(testpts,xyz,tess,tol)

% get array sizes (m points, p dimensions)
p=size(xyz,2);
[n,c]=size(testpts);
if p~=c
 error 'testpts and xyz must have the same number of columns'
end
if p<2
 error 'Points must lie in at least a 2-d space.'
end

% was the convex hull supplied?
if (nargin<3) || isempty(tess)
 tess=convhulln(xyz);
end
[nt,c]=size(tess);
if c~=p
 error 'tess array is incompatible with a dimension p space'
end

% was tol supplied?
if (nargin<4) || isempty(tol)
 tol=0;
end

% build normal vectors
switch p
 case 2
 % really simple for 2-d
  nrmls=(xyz(tess(:,1),:)-xyz(tess(:,2),:))*[0 1;-1 0];
    
  % Any degenerate edges?
  del=sqrt(sum(nrmls.^2,2));
  degenflag=(del<(max(del)*10*eps));
  if sum(degenflag)>0
   warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
   ' degenerate edges identified in the convex hull'])
   % we need to delete those degenerate normal vectors
   nrmls(degenflag,:)=[];
   nt=size(nrmls,1);
  end
  case 3
   % use vectorized cross product for 3-d
   ab=xyz(tess(:,1),:)-xyz(tess(:,2),:);
   ac=xyz(tess(:,1),:)-xyz(tess(:,3),:);
   nrmls=cross(ab,ac,2);
   degenflag=false(nt,1);
  otherwise
   % slightly more work in higher dimensions, 
   nrmls=zeros(nt,p);
   degenflag=false(nt,1);
   for i=1:nt
    % just in case of a degeneracy
    nullsp=null(xyz(tess(i,2:end),:)-repmat(xyz(tess(i,1),:),p-1,1))';
    if size(nullsp,1)>1
     degenflag(i)=true;
     nrmls(i,:)=NaN;
    else
     nrmls(i,:) = nullsp;
    end
   end
   if sum(degenflag)>0
    warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
        ' degenerate simplexes identified in the convex hull'])
    % we need to delete those degenerate normal vectors
    nrmls(degenflag,:)=[];
    nt=size(nrmls,1);
   end
end

% scale normal vectors to unit length
nrmllen=sqrt(sum(nrmls.^2,2));
nrmls=nrmls.*repmat(1./nrmllen,1,p);

% center point in the hull
center=mean(xyz,1);

% any point in the plane of each simplex in the convex hull
a=xyz(tess(~degenflag,1),:);

% ensure the normals are pointing inwards
dp=sum((repmat(center,nt,1)-a).*nrmls,2);
k=dp<0;
nrmls(k,:)=-nrmls(k,:);

% We want to test if:  dot((x - a),N) >= 0. If so for all faces of the 
% hull, then x is inside the hull. Change this to dot(x,N) >= dot(a,N)
aN=sum(nrmls.*a,2);

% test, be careful in case there are many points
in=false(n,1);

% if n is too large, the dot product may grab huge chunks of memory.
memblock=1e6;
blocks=max(1,floor(n/(memblock/nt)));
aNr=repmat(aN,1,length(1:blocks:n));
for i=1:blocks
 j=i:blocks:n;
 if size(aNr,2)~=length(j)
  aNr=repmat(aN,1,length(j));
 end
 in(j)=all((nrmls*testpts(j,:)'-aNr)>=-tol,1)';
end

