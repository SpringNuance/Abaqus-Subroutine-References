% Geometry of a planar polygon
% Based on the original function by Sommer III

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function [geom]=polygeom(x,y)

% check if inputs are same size
if ~isequal(size(x),size(y))
  error('X and Y must be the same size');
end

% number of vertices
[x,~]=shiftdim(x);
[y,~]=shiftdim(y);
[n,~]=size(x);

% temporarily shift data to mean of vertices for improved accuracy
xm = mean(x);
ym = mean(y);
x = x - xm*ones(n,1);
y = y - ym*ones(n,1);

% delta x and delta y
dx = x( [ 2:n 1 ] ) - x;
dy = y( [ 2:n 1 ] ) - y;

% summations for CW boundary integrals
A = sum( y.*dx - x.*dy )/2;
Axc = sum( 6*x.*y.*dx -3*x.*x.*dy +3*y.*dx.*dx +dx.*dx.*dy )/12;
Ayc = sum( 3*y.*y.*dx -6*x.*y.*dy -3*x.*dy.*dy -dx.*dy.*dy )/12;
P = sum( sqrt( dx.*dx +dy.*dy ) );

% check for CCW versus CW boundary
if A < 0
  A = -A;
  Axc = -Axc;
  Ayc = -Ayc;
end

% centroidal moments
xc = Axc / A;
yc = Ayc / A;

% replace mean of vertices
x_cen = xc + xm;
y_cen = yc + ym;

% return values
geom=[A x_cen y_cen P];


