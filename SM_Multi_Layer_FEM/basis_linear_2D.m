% function [value,d_value]=basis_linear_2D(x)
% M=length(x); % 2X3
% value=zeros(3,M);
% value(1,:)=ones(1,M)-x(1,:)-x(2,:); % Row Vectors
% value(2,:)=x(1,:);
% value(3,:)=x(2,:);
% d_value=zeros(2,M,3);
% v=ones(1,M);
% d_value(:,:,1)=[-v;-v];
% d_value(:,:,2)=[v;zeros(1,M)];
% d_value(:,:,3)=[zeros(1,M);v];
% end

function [value,d_value]=basis_linear_2D(x)
M = size(x,2); % number of quadrature points
value = zeros(3,M);
value(1,:) = 1 - x(1,:) - x(2,:); % φ1 = 1 - ξ - η
value(2,:) = x(1,:);              % φ2 = ξ
value(3,:) = x(2,:);              % φ3 = η

d_value = zeros(2,M,3);
v = ones(1,M);
d_value(:,:,1) = [-v; -v];       % ∇φ1 = (-1, -1)
d_value(:,:,2) = [ v;  0*v];     % ∇φ2 = (1, 0)
d_value(:,:,3) = [ 0*v;  v];     % ∇φ3 = (0, 1)
end
