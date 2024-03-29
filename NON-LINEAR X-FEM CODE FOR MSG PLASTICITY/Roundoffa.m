% Function to round a scalar, matrix or vector

% If using this code for research or industrial purposes please cite:
% E. Martínez-Pañeda, S. Natarajan, S. Bordas. 
% Gradient plasticity crack tip characterization by means of the extended
% finite element method. Computational Mechanics (2017)
% doi:10.1007/s00466-017-1375-6

function y = Roundoffa(number,decimal_places)

[INeg,JNeg]=find(number<0); % Negative numbers

if ~isempty(INeg)
 IndNeg = sub2ind(size(number),INeg,JNeg);
 Number = abs(number);
else
 Number = number;
end

decimals = 10.^decimal_places;
y1 = fix(decimals * Number + 0.5)./decimals;

if ~isempty(INeg)
 y1(IndNeg) = -y1(IndNeg);
end

y = y1;