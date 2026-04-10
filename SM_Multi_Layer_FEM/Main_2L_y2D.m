clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
schnakenberg=01;
hopf=01;
boundary=0;
out_of_phase=0;
in_phase=0;
savemat=1;

drawperframe=500;
tol = 1e-9;

%% time discretization
T = 250;
dt = 0.001;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%%
if savemat
if abs(T-250) < 1e-8
    saveStride = 1;
else
    saveStride = drawperframe;   % fallback
end
nsave = ceil(nt / saveStride);
end

%%
N = 2;
L = 20;
H = 0.1;
eta = 0.1; 

%% saving folder and name
folder='D:\MATLAB Output\SM_Multi_Layer_FEM\';

%% Define Parameters, the Source Term, Coupling Term and steady states
if schnakenberg
    if out_of_phase
folder = [folder,'out_of_phase\'];
a=0.05; b=1.2;
D_2=1;
D_1=30;    
    elseif hopf
folder = [folder,'hopf\'];
a=0.05; b=0.5;
D_2=1;
D_1=1;
    elseif boundary
folder = [folder,'boundary\'];
a=0.05; b=0.7;
D_2=1;
D_1=20;
    end        
f_2=@(u) a - u;
f_1=@(v) b - 0*v;
G_2=@(u2,u1) (u2.^2.*u1);
G_1=@(u1,u2) (u2.^2.*u1);
u20=a+b; u10=b/(a+b)^2;
sign = [+1;-1];
elseif in_phase
folder = [folder,'in_phase\'];
a=0.75; b=-1; ep=0.01;
D_2=1;
D_1=50;
f_2=@(u) a*u - ep*u.^3;
f_1=@(v) b*v - ep*v.^3;
G_2=@(u2,u1) (u1 + u2);
G_1=@(u1,u2) (u1 + u2);
u20=0; u10=0;
sign = [-1;+1];
end

%%
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
diary OFF;
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

%%
FileNameAndLocation = mfilename('fullpath');
[filepath,name,~] = fileparts(FileNameAndLocation);
ext='.m';
origFile = fullfile(filepath,[name ext]);
newbackup = fullfile(folder, sprintf('%s_%s.txt',time,name));
copyfile(origFile, newbackup);
fprintf('Backup created: %s\n', newbackup);

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

whos p1 t1

n2 = size(p2,2);
n1 = size(p1,2);

L2=[0,L,H,2*H];
L1=[0,L,0,H];

%% Steady states
fprintf(['\nInitial steady state for\n' ...
    'surface layer is: u20=%.4f\n'...
    'core layer is: u10=%.4f\n'], ...
    u20,u10);

%% Assemble Matrices
ord=3;
[S_2,M_2]=AssembleGlobalMatrices2D(p2,t2,ord);
[S_1,M_1]=AssembleGlobalMatrices2D(p1,t1,ord);

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

u2_range = zeros(nt, 1);
u1_range = zeros(nt, 1);

%%
for e=1:length(eta)
evall=round(eta(e),4);
fprintf('\nFor eta=%.4f,\n',evall);
eval=['_',num2str(evall)];

%%
if savemat
bigMatName0 = [prefix, eval, '_timeseries0.mat'];
bigMatName1 = [prefix, eval, '_timeseries1.mat'];
bigMatName2 = [prefix, eval, '_timeseries2.mat'];

mobj0 = matfile(bigMatName0,'Writable',true);
mobj1 = matfile(bigMatName1,'Writable',true);
mobj2 = matfile(bigMatName2,'Writable',true);

mobj0.tt  = zeros(nsave,1,'double');
if hopf
mobj0.u1_avg = zeros(nsave,1,'double');
mobj0.u2_avg = zeros(nsave,1,'double');
end
mobj1.u1  = zeros(n1,nsave,'single');
mobj2.u2  = zeros(n2,nsave,'single');
end

%%
pattern_detected = false;
pattern_end=5;
stopti=nt; 

%% Initial condition
rng(0); 
Perturbations2 = 0.01*randn(n2,1);
Perturbations1 = 0.01*randn(n1,1);
u2 = u20 + Perturbations2; 
u1 = u10 + Perturbations1; 

%% set up figure
vq2=griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
vq1=griddata(p1(1,:),p1(2,:),u1,xq1,yq1);

giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');
ax2 = axes('Parent',fig);   % for top layer
hold(ax2,'on');
u2_fig=surf(ax2, xq2, yq2, vq2, 'EdgeColor','none');
view(ax2,2)
xlim(ax2,[L2(1),L2(2)])
ylim(ax2,[L2(3),L2(4)])
axis(ax2,'tight')
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$u_2,~y$','Interpreter','latex')
ax2.Position=[0.1 0.51 0.74 0.4];
set(ax2,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax2.FontSize=12;
ax2.YTick=[L2(3) L2(4)];
colormap(ax2,'jet')
cb2=colorbar(ax2);
cb2.Position=[0.87 0.512 0.02 0.398];
cb2.Label.Units='normalized';
cb2.FontSize = 10;

ax1 = axes('Parent',fig);   % for bottom layer
hold(ax1,'on');
u1_fig=surf(ax1, xq1, yq1, vq1, 'EdgeColor','none');
view(ax1,2)
xlim(ax1,[L1(1),L1(2)])
ylim(ax1,[L1(3),L1(4)])
axis(ax1,'tight')
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$u_1,~y$','Interpreter','latex')
ax1.Position=[0.1 0.11 0.74 0.4];
ax1.FontSize=12;
ax1.YTick=[L1(3) L1(4)];
colormap(ax1,'jet')
cb1=colorbar(ax1);
cb1.Position=[0.87 0.11 0.02 0.398];
cb1.Label.Units='normalized';
cb1.FontSize = 10;

saveID = 0;

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

    %% --- Right-hand side assembly ---
    ord=3;
    F_2 = ReactKineInt(p2,t2,u2,f_2,ord);
    F_1 = ReactKineInt(p1,t1,u1,f_1,ord);

    ord=2;
    b2_bottom = NLBoundFluxInt(bottom_t2, top_t1, p2, p1, u2, u1, G_2, sign(1)*eta(e), ord);
    b1_top    = NLBoundFluxInt(top_t1, bottom_t2, p1, p2, u1, u2, G_1, sign(2)*eta(e), ord);

    %% Matrix system
    u2_new = (M_2 + dt/2*D_2*S_2)\ ...
            ((M_2 - dt/2*D_2*S_2)*u2 + dt*(F_2 + b2_bottom));

    u1_new = (M_1 + dt/2*D_1*S_1)\ ...
            ((M_1 - dt/2*D_1*S_1)*u1 + dt*(F_1 + b1_top));

    %% Pattern detection using norms
    err2 = norm(u2_new - u2)/norm(u2);
    err1 = norm(u1_new - u1)/norm(u1);
    maxErr = max([err2, err1]);

    %% Update and save
    u2 = u2_new;
    u1 = u1_new;

    % u2_avg(ti) = mean(u2);
    % u1_avg(ti) = mean(u1);

    u2_avg(ti) = (ones(1,n2) * M_2 * u2) / (ones(1,n2) * M_2 * ones(n2,1));
    u1_avg(ti) = (ones(1,n1) * M_1 * u1) / (ones(1,n1) * M_1 * ones(n1,1));

    % Calculate max-min for each layer
    u2_range(ti) = max(u2) - min(u2);
    u1_range(ti) = max(u1) - min(u1);
    
    %%
    if savemat
    saveID = saveID + 1;
    mobj0.tt(saveID,1) = t;
    if hopf
    mobj0.u1_avg(saveID,1) = u1_avg(ti);
    mobj0.u2_avg(saveID,1) = u2_avg(ti);    
    end
    mobj1.u1(:,saveID) = u1;
    mobj2.u2(:,saveID) = u2;
    end

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        vq2 = griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
        vq1 = griddata(p1(1,:),p1(2,:),u1,xq1,yq1);
        if showanimation
            u2_fig.ZData = vq2;
            u1_fig.ZData = vq1;
            ax2.Title.String = sprintf('t = %.4f', t);
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
%%
if savemat
mobj0.p1 = p1;  mobj0.p2 = p2;  

mobj0.L1 = L1;  mobj0.L2 = L2;  

mobj0.xq1 = xq1;  mobj0.yq1 = yq1;
mobj0.xq2 = xq2;  mobj0.yq2 = yq2;

mobj0.dt = dt;
mobj0.T = T;
mobj0.nsave = saveID;
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
xlabel('$t$','Interpreter','latex');
ylabel('Spatial average','Interpreter','latex');
legend('$\bar u_2(t)$', '$\bar u_1(t)$', 'Interpreter','latex');
title('Spatially averaged dynamics (Hopf)');
grid on;
saveas(fig_avg,[prefix,eval,'_avg.png']);
saveas(fig_avg,[prefix,eval,'_avg.fig']);

%% Plot Spatial Range (Max - Min)
u1_range_used = u1_range(1:ti);
u2_range_used = u2_range(1:ti);

fig_range = figure;
plot(t_used, u2_range_used, 'r.', 'LineWidth', 2); hold on;
plot(t_used, u1_range_used, 'b--', 'LineWidth', 2);
xlabel('$t$','Interpreter','latex');
ylabel('$\max(u) - \min(u)$', 'Interpreter', 'latex');
legend('Range $u_2$ (Top)', 'Range $u_1$ (Bottom)', 'Interpreter', 'latex');
title('Pattern Amplitude (Spatial Max - Min) over Time');
ax = gca;
ax.FontSize=12;
grid on;
saveas(fig_range, [prefix, eval, '_range.png']);
saveas(fig_range, [prefix, eval, '_range.fig']);

%% ODE
if schnakenberg
odefun = @(t,y) [
    b - eta(e)/H * y(2)^2 * y(1);
    a - y(2) + eta(e)/H * y(2)^2 * y(1)
];
elseif in_phase
odefun = @(t,y) [
    b * y(1) - ep * y(1)^3 + eta(e)/H * (y(1) + y(2));
    a * y(2) - ep * y(2)^3 - eta(e)/H * (y(1) + y(2));
];
end

y0 = [u10; u20] + 1e-3* [1; -1];
[t_ode, y_ode] = ode45(odefun, [0 t], y0);

%%
fig_ode=figure;
plot(t_ode, y_ode(:,1), 'b--', 'LineWidth', 2); hold on;
plot(t_ode, y_ode(:,2), 'r--', 'LineWidth', 2);
xlabel('$t$','Interpreter','latex');
ylabel('Concentration','Interpreter','latex');
legend('ODE $u_1$', ...
       'ODE $u_2$', ...
       'Interpreter','latex');
title('Hopf oscillations: ODE');
grid on;
saveas(fig_ode,[prefix,eval,'_ode.png']);
saveas(fig_ode,[prefix,eval,'_ode.fig']);

%%
fig_ode_pde=figure;
plot(t_used, u1_used, 'b', 'LineWidth', 2); hold on;
plot(t_ode, y_ode(:,1), 'b--', 'LineWidth', 2); hold on;
plot(t_used, u2_used, 'r', 'LineWidth', 2);
plot(t_ode, y_ode(:,2), 'r--', 'LineWidth', 2);
xlabel('$t$','Interpreter','latex');
ylabel('Concentration','Interpreter','latex');
legend('$\bar u_1$ (2D PDE)', '$u_1$ (ODE)', ...
       '$\bar u_2$ (2D PDE)', '$u_2$ (ODE)', ...
       'Interpreter','latex');
title('Hopf oscillations: PDE vs ODE');
ax = gca;
ax.FontSize=12;
grid on;
saveas(fig_ode_pde,[prefix,eval,'_ode_pde.png']);
saveas(fig_ode_pde,[prefix,eval,'_ode_pde.fig']);

end

%%
fprintf('\nDone!\n');
toc
diary OFF
