% Function to force an integer to lie between given limits by wrapping
% Based on the original function by John Burkardt

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function value=i4_wrap(ival,ilo,ihi)

jlo=min(ilo,ihi);
jhi=max(ilo,ihi);

wide=jhi-jlo+1;

if(wide==1)
 value=jlo;
else
 i=ival-jlo;
 j=wide;
 if(j==0)
  fprintf(1,'\n');
  fprintf(1,'i4_wrap - Fatal error!\n');
  fprintf(1,'Illegal divisor J = %d\n',j);
  error('i4_wrap - Fatal error!');
 end
 value1=mod(i,j);
 if (value1<0)
  value1=value1+abs(j);
 end
 value=jlo+value1;
end
