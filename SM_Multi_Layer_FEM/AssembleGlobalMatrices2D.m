function [St,Ma] = AssembleGlobalMatrices_corrected2(p,t,ord)
% Assemble global stiffness (St) and mass (Ma) matrices
% p: 2xN node coordinates
% t: 3xM element connectivity
% ord: quadrature order

N = size(p,2);   % number of nodes
M = size(t,2);   % number of elements

St = sparse(N,N);
Ma = sparse(N,N);

%% Quadrature rule and basis functions on reference element
[iw,ip] = RefElemQuad(ord);           % iw: 1xnq, ip: 2xnq
[basis_ip,grad_basis_ip] = basis_linear_2D(ip);
% basis_ip(i,:) = values of basis i at quadrature points
% grad_basis_ip(:,q,i) = gradient of basis i at quadrature point q

%% Loop over all triangles
for k = 1:M
    nodes = p(:,t(:,k));              % 2x3 coordinates of element
    C = [nodes(:,2)-nodes(:,1), nodes(:,3)-nodes(:,1)]; % Jacobian
    Cdet = det(C);
    CinvT = inv(C)';                  % (C^-1)^T

    node_ids = t(:,k);                % 3 local node indices
    Local_M = zeros(3,3);
    Local_S = zeros(3,3);

    % Loop over basis functions
    for i = 1:3
        for j = 1:3
            % Mass matrix entry
            Local_M(i,j) = (basis_ip(i,:) .* basis_ip(j,:)) * iw' * Cdet;

            % Stiffness matrix entry (integral of grad phi_i · grad phi_j)
            gi = CinvT * squeeze(grad_basis_ip(:,:,i)); % 2xnq
            gj = CinvT * squeeze(grad_basis_ip(:,:,j)); % 2xnq
            Local_S(i,j) = sum(sum(gi .* gj,1) .* iw) * Cdet;
        end
    end

    % Assemble
    St(node_ids,node_ids) = St(node_ids,node_ids) + Local_S;
    Ma(node_ids,node_ids) = Ma(node_ids,node_ids) + Local_M;
end
end
