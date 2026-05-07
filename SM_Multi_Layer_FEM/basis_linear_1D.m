function [basis, grad_basis] = basis_linear_1D(x)
% Evaluates linear basis functions and their gradients on [0,1]
x = x(:)';  % Ensure row vector for consistent output
n = length(x);

basis = [1 - x; x];        % phi1 = 1 - x, phi2 = x

grad_basis = repmat([-1; 1], 1, n);  % dphi1 = -1, dphi2 = 1
end
