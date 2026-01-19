function b = NLBoundFluxInt(edges1, edges2, p1, p2, u1, u2, G, eta, ord)
% Computes 1D boundary flux contribution
% Number of nodes
N = size(p1, 2);
b = zeros(N, 1);

% Quadrature rule and basis functions
[iw_bc, ip_bc] = RefEdgeQuad(ord);
[basis_ip_bc,~] = basis_linear_1D(ip_bc);
% size(edges1,1)
% size(edges1,2)
% Loop over edges
for e = 1:size(edges1, 2)
    nodes1 = edges1(:, e); %whos nodes1 edges1
    nodes2 = edges2(:, e); %whos nodes2 edges2
    x1 = p1(:, nodes1);
    x2 = p2(:, nodes2);
    % x1-x2
    % size(nodes1)
    % none
    % isequal(x1,x2)
    edge_vec = x1(:,2) - x1(:,1);
    edge_length = norm(edge_vec);

    for q = 1:numel(ip_bc)
        phi = basis_ip_bc(:, q);        % 2x1 vector at this quad point
%         size(phi)
% size(u1(nodes1))
% size(u2(nodes2))
        % Interpolate u1 and u2 at quadrature point
        u1_q = phi' * u1(nodes1);
        u2_q = phi' * u2(nodes2);

        % Evaluate nonlinear flux
        g_val = eta * G(u1_q, u2_q);
        
if any(isnan(u1_q)) || any(isnan(u2_q))
    error('NaN detected in u_ip values at edge %d (quadrature point %d)', e, q);
end

if any(isnan(g_val))
    disp('Inputs to g:');
    disp([u1_q; u2_q]);
    error('G returned NaN at edge %d', e);
end
if any(isinf(u1_q)) || any(isinf(u2_q))
    error('Inf detected in u_ip values at edge %d', e);
end

if any(isinf(g_val))
    disp('Inputs to g:');
    disp([u1_q; u2_q]);
    error('G returned Inf at edge %d', e);
end

        % Accumulate into global RHS vector
        b(nodes1) = b(nodes1) + iw_bc(q) * edge_length * g_val * phi;
    end
end
% [u1_q u2_q]
end
