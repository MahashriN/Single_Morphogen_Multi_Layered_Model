clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
tol = 1e-9;

T = 200;
dt = 0.001;

N = 2;
L = 20;
H = 0.1;
eta = 0.25; 

%% saving folder and name
folder='D:\MATLAB Output\SM_Multi_Layer_FEM\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
diary OFF;
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

n2 = size(p2,2);
n1 = size(p1,2);

L2=[0,L,H,2*H];
L1=[0,L,0,H];

%% Parameters
a=0.05; b=0.7;
D_2=1;
D_1=20;

%% Define the Source Term and Coupling Term
f_2=@(u) a - u;
f_1=@(v) b - 0*v;

G_2=@(u2,u1) u2.^2.*u1;
G_1=@(u1,u2) u2.^2.*u1;

%% Steady states
u20=a+b; u10=b/(a+b)^2;

fprintf(['\nInitial steady state for\n' ...
    'surface layer is: u20=%.4f\n'...
    'core layer is: u10=%.4f\n'], ...
    u20,u10);

%% time discretization
% T=1000;
% dt=0.01; %1e-4;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%% Assemble Matrices
ord=3;tic
[S_2,M_2]=AssembleGlobalMatrices2D(p2,t2,ord);
[S_1,M_1]=AssembleGlobalMatrices2D(p1,t1,ord);
M_big = blkdiag(M_2,M_1);
S_big = blkdiag(D_2*S_2, D_1*S_1);
toc

%% Grids
nx = 200; ny = 50;

[xq2,yq2] = meshgrid(linspace(L2(1),L2(2),nx), ...
                     linspace(L2(3),L2(4),ny));

[xq1,yq1] = meshgrid(linspace(L1(1),L1(2),nx), ...
                     linspace(L1(3),L1(4),ny));

%% Preallocate storage
max_steps = nt;
tt = zeros(max_steps,1); 

u2_avg = zeros(nt,1);
u1_avg = zeros(nt,1);

%%
for e=1:length(eta)
evall=round(eta(e),4);
fprintf('\nFor eta=%.4f,\n',evall);
eval=['_',num2str(evall)];

pattern_detected = false;
pattern_end=5;
stopti=nt; 

%% Initial condition
rng(0); 
Perturbations2 = 1+0.1*randn(n2,1);
Perturbations1 = 1+0.1*randn(n1,1);
u2=u20 * Perturbations2;
u1=u10 * Perturbations1;

%% set up figure
vq2=griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
vq1=griddata(p1(1,:),p1(2,:),u1,xq1,yq1);

giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');
ax2 = axes('Parent',fig);   % for u (top layer)
hold(ax2,'on');
u2_fig=surf(ax2, xq2, yq2, vq2, 'EdgeColor','none');
view(ax2,2)
xlim(ax2,[L2(1),L2(2)])
ylim(ax2,[L2(3),L2(4)])
axis(ax2,'tight')
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$u_2,~y$','Interpreter','latex')
ax2.Position=[0.1 0.51 0.67 0.4];
set(ax2,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax2.FontSize=12;
ax2.YTick=[L2(3) L2(4)];
colormap(ax2,'jet')
cb2=colorbar(ax2);
cb2.Position=[0.9 0.11 0.015 0.8];
cb2.Label.String = '$u_2$';
cb2.Label.Rotation=0;
cb2.Label.Interpreter='latex';
cb2.Label.Units='normalized';
cb2.Label.Position = [0.6 -0.04 0];
cb2.FontSize = 10;

ax1 = axes('Parent',fig);   % for v (bottom layer)
hold(ax1,'on');
u1_fig=surf(ax1, xq1, yq1, vq1, 'EdgeColor','none');
view(ax1,2)
xlim(ax1,[L1(1),L1(2)])
ylim(ax1,[L1(3),L1(4)])
axis(ax1,'tight')
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$u_1,~y$','Interpreter','latex')
ax1.Position=[0.1 0.11 0.67 0.4];
ax1.FontSize=12;
ax1.YTick=[L1(3) L1(4)];
colormap(ax1,'jet')
cb1=colorbar(ax1);
cb1.Position=[0.81 0.11 0.015 0.8];
cb1.Label.String = '$u_1$';
cb1.Label.Rotation=0;
cb1.Label.Interpreter='latex';
cb1.Label.Units='normalized';
cb1.Label.Position = [0.6 -0.04 0];
cb1.FontSize = 10;

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

%% --- Right-hand side assembly ---
ord=3;
F_2 = ReactKineInt(p2,t2,u2,f_2,ord);
F_1 = ReactKineInt(p1,t1,u1,f_1,ord);

ord=2;
b2_bottom = NLBoundFluxInt(bottom_t2, top_t1, p2, p1, u2, u1, G_2, eta(e), ord);
b1_top    = NLBoundFluxInt(top_t1, bottom_t2, p1, p2, u1, u2, G_1, -eta(e), ord);

if any(isinf(b2_bottom(:))) || any(isinf(b1_top(:)))
    error('Inf at the boundary conditions  at step %d (time = %.5f)', ti, t);
end

%% Matrix system
u2_new = (M_2 + dt/2*D_2*S_2)\ ...
         ((M_2 - dt/2*D_2*S_2)*u2 + dt*(F_2 + b2_bottom));

u1_new = (M_1 + dt/2*D_1*S_1)\ ...
         ((M_1 - dt/2*D_1*S_1)*u1 + dt*(F_1 + b1_top));

%%
U_big = [u2; u1];
F_big = [F_2; F_1];
B_big = [b2_bottom; b1_top];
A = M_big + 0.5*dt*S_big;
B = (M_big - 0.5*dt*S_big)*U_big + dt*(F_big + B_big);

%% Stability and NaN check
if any(isnan(A(:))) || any(isnan(B(:)))
 [any(isnan(A(:)))   %1
any(isnan(B(:)))    %2
any(isnan(u2(:)))   %3
any(isnan(u1(:)))   %4
any(isnan(B_big(:)))    %5
any(isnan(F_big(:)))    %6
any(isnan(b2_bottom(:)))    %7
any(isnan(b1_top(:)))]    %8
   error('Matrix A or RHS B contains NaN at step %d (time = %.5f)', ti, t);
end
if condest(A) > 1e12
    warning('Matrix A is ill-conditioned at step %d (time = %.5f), condest = %.2e', ti, t, condest(A));
end
    %% Pattern detection using norms
    err2 = norm(u2_new - u2)/norm(u2);
    err1 = norm(u1_new - u1)/norm(u1);
    maxErr = max([err2, err1]);

%% Update and save
u2 = u2_new;
u1 = u1_new;

u2_avg(ti) = mean(u2);
u1_avg(ti) = mean(u1);

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        vq2 = griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
        vq1 = griddata(p1(1,:),p1(2,:),u1,xq1,yq1);
        if showanimation
            u2_fig.ZData = vq2;
            u1_fig.ZData = vq1;
ax2.Title.String = sprintf('t = %.4f,  \\eta = %.3g', t, evall);
ax2.Title.Interpreter = 'tex';
zline = max([vq1(:); vq2(:)]) * ones(1,2);
plot3(ax2, [L1(1) L1(2)], [L2(3) L2(3)], zline, 'k','LineWidth',1);
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
t_used = tt(1:ti);
u1_used = u1_avg(1:ti);
u2_used = u2_avg(1:ti);

fig_avg=figure;
plot(t_used, u2_used, 'r', 'LineWidth', 2); hold on;
plot(t_used, u1_used, 'b', 'LineWidth', 2);
xlabel('Time');
ylabel('Spatial average');
legend('$\bar u_2(t)$', '$\bar u_1(t)$', 'Interpreter','latex');
title('Spatially averaged dynamics (Hopf)');
grid on;
saveas(fig_avg,[prefix,eval,'_avg.png']);
saveas(fig_avg,[prefix,eval,'_avg.fig']);

%%
odefun = @(t,y) [
    b - eta * y(2)^2 * y(1);
    a - y(2) + eta * y(2)^2 * y(1)
];

y0 = [u10; u20];  % same steady state
[t_ode, y_ode] = ode45(odefun, [0 t], y0);

fig_ode=figure;
plot(t_used, u1_used, 'b', 'LineWidth', 2); hold on;
plot(t_ode, y_ode(:,1), 'b--', 'LineWidth', 2);
plot(t_used, u2_used, 'r', 'LineWidth', 2);
plot(t_ode, y_ode(:,2), 'r--', 'LineWidth', 2);
xlabel('Time');
ylabel('Concentration');
legend('PDE $\bar u_1$', 'ODE $u_1$', ...
       'PDE $\bar u_2$', 'ODE $u_2$', ...
       'Interpreter','latex');
title('Hopf oscillations: PDE vs ODE');
grid on;
saveas(fig_ode,[prefix,eval,'_ode.png']);
saveas(fig_ode,[prefix,eval,'_ode.fig']);

end

%%
fprintf('\nDone!\n');
toc
diary OFF