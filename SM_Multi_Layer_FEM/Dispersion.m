folder='D:\\MATLAB Output\\SM_Multi_Layer_FEM\\Turing_wave\\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

d=0;
p1 = -3-d;p2 = -5-d;p3 = -1-d;
p4 = 8; p5 = 10; p6 = 4; p7=-2; 
delta=0.1;
eta = 0.05; %0.1;
H=0.05;

D1 = delta; D2 = 1; D3 = delta;

eta_hat=eta/H;

J = diag([p1, p2, p3]) + eta_hat*[p4, p5, 0; -p4, -p5+p6, p7; 0, -p6, -p7];

D = diag([D1, D2, D3]);

eigfun = @(mu) eig(J - mu*D);
mus = linspace(0,200,1500);

close all;
fig=figure;
hold on
for i = 1:length(mus)
    ev = eigfun(mus(i));
    [max_val, idx] = max(real(ev));
    alpha_max(i) = max(real(ev));
    omega_max(i) = imag(ev(idx));
    % filter values below -2
    ev = ev(real(ev) >= -2);
    if isempty(ev), continue; end
    used = false(size(ev));
    for a = 1:length(ev)
        if used(a), continue; end
        va = ev(a);
        if abs(imag(va)) < 1e-8
            % real eigenvalue -> black
            plot(mus(i), real(va), 'ko', 'MarkerSize', 4);
            used(a) = true;
        else
            % complex (part of a conjugate pair) -> blue
            % find its conjugate pair
            conjIdx = find(~used & abs(ev - conj(va)) < 1e-6, 1);
            rl(i)=real(va);
            if isempty(conjIdx)
                % stray complex (plot as blue real part)
                plot(mus(i), real(va), 'o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 4, 'LineWidth', 1.5);
                used(a) = true;
            else
                % plot both members of the conjugate pair as blue (same x)
                plot(mus(i), real(va), 'o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 4, 'LineWidth', 1.5);
                plot(mus(i), real(ev(conjIdx)), 'o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 4, 'LineWidth', 1.5);
                used(a) = true;
                used(conjIdx) = true;
            end
        end
    end
end
yline(0,'k--','LineWidth',1);
box on
xlabel('$k^2$','Interpreter','latex');
ylabel('$\lambda$','Interpreter','latex');
h_real = plot(nan, nan, 'ko', 'MarkerSize', 4, 'LineWidth', 1.5);
h_complex = plot(nan, nan, 'o', 'Color', [0 0.4470 0.7410], 'MarkerSize', 4, 'LineWidth', 1.5);
legend([h_real, h_complex], {'Real', 'Complex'}, 'Interpreter', 'latex');
ax=gca;
ax.FontSize=16;
saveas(fig,[prefix,sprintf('turing_wave_DC_H_%.2f_eta_%.2f.png',H,eta)])
saveas(fig,[prefix,sprintf('turing_wave_DC_H_%.2f_eta_%.2f.fig',H,eta)])
