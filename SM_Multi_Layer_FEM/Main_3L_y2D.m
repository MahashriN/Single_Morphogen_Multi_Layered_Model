diary OFF;
clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
pattern_detected = false;
tol=1e-8;

T=500;
dt=0.001;
time_snaps=50001;

N=3;
L=20;
H=0.1;
eta_23 = 0.1; 
eta_12 = 0.1; 

%% saving folder and name
folder='D:\MATLAB Output\SM_Multi_Layer_FEM\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

%%
FileNameAndLocation = mfilename('fullpath');   % current script full path without extension
[filepath,name,~] = fileparts(FileNameAndLocation);
ext='.m';
% Original script full filename
origFile = fullfile(filepath,[name ext]);

% Backup filename (adds "backup" + version + .txt)
newbackup = fullfile(folder, sprintf('%s_%s.txt',time,name));

% Check if backup already exists
A = exist(newbackup,'file');
if (A ~= 0)
    warning('Backup already exists for the current version')
else
    % Create backup by copying the current .m file
    copyfile(origFile, newbackup);
    fprintf('Backup created: %s\n', newbackup);
end

%% Load Mesh Data: p-nodes; t-elements;
folderpath = [filepath,'\N=',num2str(N),'\L=',num2str(L),'\H=',num2str(H)]; 
filePattern = fullfile(folderpath, '*.mat');
matFiles = dir(filePattern);

numFiles = length(matFiles);
for k = 1:numFiles
    currentFileName = fullfile(matFiles(k).folder, matFiles(k).name);
    fprintf('Processing file: %s \n', currentFileName);
    load(currentFileName)
end

n3 = size(p3,2);
n2 = size(p2,2);
n1 = size(p1,2);

L3=[0,L,2*H,3*H];
L2=[0,L,H,2*H];
L1=[0,L,0,H];

%% Parameters
D_3 = 0.09; D_2 = 1; D_1 = 0.09;
% par=[-2.53182 -5.50353 -0.278245 30*6.41225 30*7.9846 30*5.8195 30*-1.06791];
par=[-2.53 -5.50 -0.28 192.37 239.53 174.59 -32.04];

%% Define the Source Term and Coupling Term
G_23 = @(u2,u3) par(6)*u2 + par(7)*u3;
G_32 = @(u3,u2) par(6)*u2 + par(7)*u3;

G_12 = @(u1,u2) par(4)*u1 + par(5)*u2;
G_21 = @(u2,u1) par(4)*u1 + par(5)*u2;

f_3 = @(u) par(3)*u - u.^3;
f_2 = @(u) par(2)*u - u.^3;
f_1 = @(u) par(1)*u - u.^3;

%% Steady states
u30 = 1;
u20 = 1;
u10 = 1;

fprintf(['\nInitial steady state for\n' ...
    'layer-3 is: u30=%.4f\n'...
    'layer-2 is: u20=%.4f\n'...
    'layer-1 is: u10=%.4f\n\n'], ...
    u30,u20,u10);

%% time discretization
% T=1000;
% dt=0.01; %1e-4;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%% Assemble Matrices
ord=3;tic
[S_3,M_3] = AssembleGlobalMatrices2D(p3,t3,ord);
[S_2,M_2] = AssembleGlobalMatrices2D(p2,t2,ord);
[S_1,M_1] = AssembleGlobalMatrices2D(p1,t1,ord);
toc

%%
nx = 200; ny = 50;   % thin layers -> small ny

[xq3,yq3] = meshgrid(linspace(L3(1),L3(2),nx), ...
                     linspace(L3(3),L3(4),ny));

[xq2,yq2] = meshgrid(linspace(L2(1),L2(2),nx), ...
                     linspace(L2(3),L2(4),ny));

[xq1,yq1] = meshgrid(linspace(L1(1),L1(2),nx), ...
                     linspace(L1(3),L1(4),ny));

%% Preallocate storage
max_steps = nt;
tt = zeros(max_steps,1); 

% u3_avg = zeros(nt,1);
% u2_avg = zeros(nt,1);
% u1_avg = zeros(nt,1);

%%
for e=1:length(eta_12)
eval12 = round(eta_12(e),4);
eval23 = round(eta_23(e),4);
fprintf('\nFor eta_12 = %.4f, eta_23 = %.4f\n', eval12, eval23);
eval = ['_', num2str(eval12), '_', num2str(eval23)];

pattern_detected = false;
pattern_end=5;
stopti=nt; 

%% Initial condition
rng(0); 
Perturbations3 = 0.1*randn(n3,1);
Perturbations2 = 0.1*randn(n2,1);
Perturbations1 = 0.1*randn(n1,1);
u3 = u30 * Perturbations3;
u2 = u20 * Perturbations2;
u1 = u10 * Perturbations1;

%%
vq3=griddata(p3(1,:),p3(2,:),u3,xq3,yq3);
vq2=griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
vq1=griddata(p1(1,:),p1(2,:),u1,xq1,yq1);

giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');
fig.Position = [320 150 705 500];
ax3 = axes('Parent',fig);   % for layer-2 (middle)
hold(ax3,'on');
u3_fig=surf(ax3, xq3, yq3, vq3, 'EdgeColor','none');
view(ax3,2)
xlim(ax3,[L3(1),L3(2)])
ylim(ax3,[L3(3),L3(4)])
axis(ax3,'tight')
xlabel(ax3,'$x$','Interpreter','latex')
ylabel(ax3,'$u_3,~y$','Interpreter','latex')
ax3.Position=[0.1 0.65 0.66 0.27];
set(ax3,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax3.FontSize=12;
ax3.YTick=[L3(3) L3(4)];
colormap(ax3,'jet')
cb3=colorbar(ax3);
cb3.Position=[0.92 0.11 0.015 0.81];
cb3.Label.String = '$u_3$';
cb3.Label.Rotation=0;
cb3.Label.Interpreter='latex';
cb3.Label.Units='normalized';
cb3.Label.Position = [0.6 -0.04 0];
cb3.FontSize = 10;

ax2 = axes('Parent',fig);   % for layer-2 (middle)
hold(ax2,'on');
u2_fig=surf(ax2, xq2, yq2, vq2, 'EdgeColor','none');
view(ax2,2)
xlim(ax2,[L2(1),L2(2)])
ylim(ax2,[L2(3),L2(4)])
axis(ax2,'tight')
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$u_2,~y$','Interpreter','latex')
ax2.Position=[0.1 0.38 0.66 0.27];
set(ax2,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax2.FontSize=12;
ax2.YTick=[L2(3) L2(4)];
colormap(ax2,'jet')
cb2=colorbar(ax2);
cb2.Position=[0.85 0.11 0.015 0.81];
cb2.Label.String = '$u_2$';
cb2.Label.Rotation=0;
cb2.Label.Interpreter='latex';
cb2.Label.Units='normalized';
cb2.Label.Position = [0.6 -0.04 0];
cb2.FontSize = 10;

ax1 = axes('Parent',fig);   % for layer-1 (bottom)
hold(ax1,'on');
u1_fig=surf(ax1, xq1, yq1, vq1, 'EdgeColor','none');
view(ax1,2)
xlim(ax1,[L1(1),L1(2)])
ylim(ax1,[L1(3),L1(4)])
axis(ax1,'tight')
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$u_1,~y$','Interpreter','latex')
ax1.Position=[0.1 0.11 0.66 0.27];
ax1.FontSize=12;
ax1.YTick=[L1(3) L1(4)];
colormap(ax1,'jet')
cb1=colorbar(ax1);
cb1.Position=[0.78 0.11 0.015 0.81];
cb1.Label.String = '$u_1$';
cb1.Label.Rotation=0;
cb1.Label.Interpreter='latex';
cb1.Label.Units='normalized';
cb1.Label.Position = [0.6 -0.04 0];
cb1.FontSize = 10;

%% Simulation iteration
for ti = 60152:nt
    t = dt * (ti - 1);
    tt(ti) = t;

%% --- Right-hand side assembly ---
ord=3;
F_3 = ReactKineInt(p3,t3,u3,f_3,ord);
F_2 = ReactKineInt(p2,t2,u2,f_2,ord);
F_1 = ReactKineInt(p1,t1,u1,f_1,ord);

ord=2;
b3_bottom = NLBoundFluxInt(bottom_t3, top_t2, ...
                          p3, p2, u3, u2, G_32, -eta_23(e), ord);
b2_top    = NLBoundFluxInt(top_t2, bottom_t3, ...
                          p2, p3, u2, u3, G_23,  eta_23(e), ord);
b2_bottom = NLBoundFluxInt(bottom_t2, top_t1, ...
                          p2, p1, u2, u1, G_21, -eta_12(e), ord);
b1_top    = NLBoundFluxInt(top_t1, bottom_t2, ...
                          p1, p2, u1, u2, G_12,  eta_12(e), ord);

B_3 = b3_bottom;
B_2 = b2_bottom + b2_top;
B_1 = b1_top;

% if mod(ti, drawperframe) == 1
% fprintf('\nFor t=%.2f\n',t)
% fprintf('12 balance = %.3e\n', sum(b1_top) + sum(b2_bottom));
% fprintf('23 balance = %.3e\n', sum(b2_top) + sum(b3_bottom));
% end

%% Matrix system
u3_new = (M_3 + dt/2*D_3*S_3) \ ...
         ((M_3 - dt/2*D_3*S_3)*u3 + dt*(F_3 + B_3));

u2_new = (M_2 + dt/2*D_2*S_2) \ ...
         ((M_2 - dt/2*D_2*S_2)*u2 + dt*(F_2 + B_2));

u1_new = (M_1 + dt/2*D_1*S_1) \ ...
         ((M_1 - dt/2*D_1*S_1)*u1 + dt*(F_1 + B_1));

%% Pattern detection using norms
    err3 = norm(u3_new - u3)/norm(u3);
    err2 = norm(u2_new - u2)/norm(u2);
    err1 = norm(u1_new - u1)/norm(u1);
    maxErr = max([err3, err2, err1]);

%% Update and save
u3 = u3_new;
u2 = u2_new;
u1 = u1_new;


% u3_avg(ti) = mean(u3);
% u2_avg(ti) = mean(u2);
% u1_avg(ti) = mean(u1);

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        vq3 = griddata(p3(1,:),p3(2,:),u3,xq3,yq3);
        vq2 = griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
        vq1 = griddata(p1(1,:),p1(2,:),u1,xq1,yq1);
% ---- interface lines (always visible) ----
zline3 = max(vq3(:)) * ones(1,2);
zline2 = max(vq2(:)) * ones(1,2);

% Interface between layer 3 and 2 (y = 2*dy)
plot3(ax3, [L3(1) L3(2)], [L3(3) L3(3)], ...
      zline3, 'k', 'LineWidth', 1.5);

% Interface between layer 2 and 1 (y = dy)
plot3(ax2, [L2(1) L2(2)], [L2(3) L2(3)], ...
      zline2, 'k', 'LineWidth', 1.5);

        if showanimation
            u3_fig.ZData = vq3;
            u2_fig.ZData = vq2;
            u1_fig.ZData = vq1;
ax3.Title.String = sprintf('t = %.4f,  \\eta_{12} = %.3g,  \\eta_{23} = %.3g', ...
                            t, eval12, eval23);
ax3.Title.Interpreter = 'tex';
            drawnow;
        end
        if makegif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ti == 1
                imwrite(imind, cm, giffile, 'gif', 'Loopcount', inf);
            else
                imwrite(imind, cm, giffile, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
            end
        end
        
    end
    if mod(ti, time_snaps) == 0
        fname = sprintf('%s%s_t_%06.1f.png', prefix, eval, t);
        saveas(fig, fname);
    end
    if ~pattern_detected && maxErr < tol
        fprintf('\nPattern fully developed at t = %.5f\n', t);
        pattern_detected = true;
        stopti = ti + round(pattern_end / dt);
    end

    if ti == stopti
        break;
    end
end
%% saving final pattern
saveas(fig,[prefix,eval,'_final.png']);
saveas(fig,[prefix,eval,'_final.fig']);

%%
% t_used = tt(1:ti);
% u1_used = u1_avg(1:ti);
% u2_used = u2_avg(1:ti);
% u3_used = u3_avg(1:ti);
% 
% fig_avg=figure;
% plot(t_used, u3_used, 'g', 'LineWidth', 2); hold on;
% plot(t_used, u2_used, 'r', 'LineWidth', 2); hold on;
% plot(t_used, u1_used, 'b', 'LineWidth', 2); hold on;
% xlabel('Time');
% ylabel('Spatial average');
% legend('$\bar u_3(t)$', '$\bar u_2(t)$', '$\bar u_1(t)$', 'Interpreter','latex');
% title('Spatially averaged dynamics (Hopf)');
% grid on;
% saveas(fig_avg,[prefix,eval,'_avg.png']);
% saveas(fig_avg,[prefix,eval,'_avg.fig']);

%%
% odefun = @(t,y) [
%     par(1)*y(1) - y(1)^3 + eta_12 * (par(4)*y(1) + par(5)*y(2));
%     par(2)*y(2) - y(2)^3 - eta_12 * (par(4)*y(1) + par(5)*y(2)) ...
%                          + eta_23 * (par(6)*y(2) + par(7)*y(3));
%     par(3)*y(3) - y(3)^3 - eta_23 * (par(6)*y(2) + par(7)*y(3))
% ];
% 
% y0 = [0.01; 0.01; 0.01];   % small perturbation
% 
% [t_ode, y_ode] = ode45(odefun, [0 T], y0);
% 
% fig_ode = figure;
% 
% plot(t_used, u1_used, 'b', 'LineWidth', 2); hold on;
% plot(t_ode, y_ode(:,1), 'b--', 'LineWidth', 2);
% 
% plot(t_used, u2_used, 'r', 'LineWidth', 2); hold on;
% plot(t_ode, y_ode(:,2), 'r--', 'LineWidth', 2);
% 
% plot(t_used, u3_used, 'g', 'LineWidth', 2); hold on;
% plot(t_ode, y_ode(:,3), 'g--', 'LineWidth', 2);
% 
% xlabel('Time');
% ylabel('Concentration');
% legend('PDE $\bar u_1$', 'ODE $u_1$', ...
%        'PDE $\bar u_2$', 'ODE $u_2$', ...
%        'PDE $\bar u_3$', 'ODE $u_3$', ...
%        'Interpreter','latex');
% 
% title('Temporal dynamics: PDE vs ODE (3-layer system)');
% grid on;
% 
% saveas(fig_ode,[prefix,eval,'_ode3.png']);
% saveas(fig_ode,[prefix,eval,'_ode3.fig']);

end

%%
fprintf('\nDone!\n');
toc
diary OFF