function F = ReactKineInt(p,t,u,f,ord)
N=length(p); % No of nodes
M=length(t); % No of elements
F=zeros(N,1); 

[iw,ip]=RefElemQuad(ord);
[basis_ip,~]=basis_linear_2D(ip);

for k=1:M
    nodes=p(:,t(1:3,k));
    C=[nodes(:,2)-nodes(:,1), nodes(:,3)-nodes(:,1)];
    Cdet=C(1,1)*C(2,2)-C(1,2)*C(2,1);
    
    node_ids = t(1:3,k);
    u_elem = u(node_ids); % local nodal values of u
    
    u_at_ip = u_elem' * basis_ip; % 1 x num_quad_pts
    f_nl = f(u_at_ip(:))';  % 1 x num_quad_pts
if any(isnan(u_at_ip))
    error('NaN detected in u_ip values at element %d', k);
end

if any(isnan(f_nl))
    disp('Inputs to f:');
    disp([u_at_ip]);
    error('f returned NaN at element %d', k);
end
if any(isinf(u_at_ip))
    error('Inf detected in u_ip values at element %d', k);
end

if any(isinf(f_nl))
    disp('Inputs to f:');
    disp([u_at_ip]);
    error('f returned Inf at element %d', k);
end

    Local_load = zeros(3,1);
    for i = 1:3
        Local_load(i) = (f_nl .* basis_ip(i,:)) * iw' * Cdet;
    end
    F(node_ids) = F(node_ids) + Local_load;
end
end
