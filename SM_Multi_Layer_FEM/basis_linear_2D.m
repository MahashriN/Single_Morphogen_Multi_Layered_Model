function [value,d_value]=basis_linear_2D(x)
M = size(x,2); % number of quadrature points
value = zeros(3,M);
value(1,:) = 1 - x(1,:) - x(2,:); % phi1 = 1 - a - b
value(2,:) = x(1,:);              % phi2 = a
value(3,:) = x(2,:);              % phi3 = b

d_value = zeros(2,M,3);
v = ones(1,M);
d_value(:,:,1) = [-v; -v];       % dphi1 = (-1, -1)
d_value(:,:,2) = [ v;  0*v];     % dphi2 = (1, 0)
d_value(:,:,3) = [ 0*v;  v];     % dphi3 = (0, 1)
end
