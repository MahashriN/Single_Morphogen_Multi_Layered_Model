clear all; close all; clc;
tic

loadGIF = 0;

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
pattern_detected = false;
tol=1e-8;

%% time discretization
T=250;
dt=0.001;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%%
if abs(T-250) < 1e-8
    saveStride = 1;
else
    saveStride = drawperframe; 
end
nsave = ceil(nt / saveStride);
j=int32(linspace(0,nFrame,4));
jj = [j(1):j(1)+10 j(2):j(2)+10 j(3)-10:j(3) j(4)-10:j(4)];
jjj = drawperframe*(jj)+1;

%%
N=3;
L=20;
H=0.01;
eta_23 = 1; 
eta_12 = 1; 

%% saving folder and name
folder='D:\MATLAB Output\SM_Multi_Layer_FEM\Turing_wave\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
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
n3 = size(p3,2);
n2 = size(p2,2);
n1 = size(p1,2);

L3=[0,L,2*H,3*H];
L2=[0,L,H,2*H];
L1=[0,L,0,H];

%% Parameters
delta = 0.1;
D_3 = delta; D_2 = 1; D_1 = delta; ep=1;
par=[-3 -5 -1 8 10 4 -2];

%% Define the Source Term and Coupling Term
G_23 = @(u2,u3) par(6)*u2 + par(7)*u3;
G_32 = @(u3,u2) par(6)*u2 + par(7)*u3;

G_12 = @(u1,u2) par(4)*u1 + par(5)*u2;
G_21 = @(u2,u1) par(4)*u1 + par(5)*u2;

f_3 = @(u) par(3)*u - ep*u.^3;
f_2 = @(u) par(2)*u - ep*u.^3;
f_1 = @(u) par(1)*u - ep*u.^3;

%% Steady states
u30 = 1;
u20 = 1;
u10 = 1;

fprintf(['\nInitial steady state for\n' ...
    'layer-3 is: u30=%.4f\n'...
    'layer-2 is: u20=%.4f\n'...
    'layer-1 is: u10=%.4f\n\n'], ...
    u30,u20,u10);

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

%% --- Updated kymograph storage (Packed) ---
nKymo = floor(nt / drawperframe) + 1;
% Space-Time (x vs t) - Averaging over y
ky_u1_p = zeros(nx, nKymo); 
ky_u2_p = zeros(nx, nKymo);
ky_u3_p = zeros(nx, nKymo);

% Space-Time (y vs t) - Averaging over x
kx_u1_p = zeros(ny, nKymo); 
kx_u2_p = zeros(ny, nKymo);
kx_u3_p = zeros(ny, nKymo);

t_kymo = zeros(1, nKymo); % To store the specific time points
k_idx = 0; % Counter for frames

%% Preallocate storage
max_steps = nt;

%%
for e=1:length(eta_12)
eval12 = round(eta_12(e),4);
eval23 = round(eta_23(e),4);
fprintf('\nFor eta_12 = %.4f, eta_23 = %.4f\n', eval12, eval23);
eval = ['_', num2str(eval12), '_', num2str(eval23)];

%%
bigMatName0 = [prefix, eval, '_timeseries0.mat'];
bigMatName1 = [prefix, eval, '_timeseries1.mat'];
bigMatName2 = [prefix, eval, '_timeseries2.mat'];
bigMatName3 = [prefix, eval, '_timeseries3.mat'];

mobj0 = matfile(bigMatName0,'Writable',true);
mobj1 = matfile(bigMatName1,'Writable',true);
mobj2 = matfile(bigMatName2,'Writable',true);
mobj3 = matfile(bigMatName3,'Writable',true);

mobj0.tt  = zeros(nsave,1,'double');

mobj1.u11  = zeros(n1,nsave,'single');
mobj2.u22  = zeros(n2,nsave,'single');
mobj3.u33  = zeros(n3,nsave,'single');

%%
pattern_detected = false;
pattern_end=5;
stopti=nt; 

%% Initial condition
rng(0); 
Perturbations3 = 0.01*randn(n3,1);
Perturbations2 = 0.01*randn(n2,1);
Perturbations1 = 0.01*randn(n1,1);
u3 = u30 * Perturbations3;
u2 = u20 * Perturbations2;
u1 = u10 * Perturbations1;

%%
vq3=griddata(p3(1,:),p3(2,:),u3,xq3,yq3);
vq2=griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
vq1=griddata(p1(1,:),p1(2,:),u1,xq1,yq1);

%% pre-figure settings
gmin = -Inf; 
gmax = inf; 

giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');
fig.Position = [43 100 1000 800];

ax3 = axes('Parent',fig);   % for layer-3 (top)
hold(ax3,'on');
u3_fig=surf(ax3, xq3, yq3, vq3, 'EdgeColor','none');
view(ax3,2)
xlim(ax3,[L3(1),L3(2)])
ylim(ax3,[L3(3),L3(4)])
clim(ax3, [gmin, gmax]);
axis(ax3,'tight')
xlabel(ax3,'$x$','Interpreter','latex')
ylabel(ax3,'$u_3,~y$','Interpreter','latex')
ax3.Position=[0.1 0.65 0.74 0.27];
set(ax3,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax3.FontSize=26;
ax3.YTick=[L3(3) L3(4)];
colormap(ax3,'jet')
ax3.YLabel.Units='normalized';
ax3.YLabel.Position=[-0.07 0.5 0];

ax2 = axes('Parent',fig);   % for layer-2 (middle)
hold(ax2,'on');
u2_fig=surf(ax2, xq2, yq2, vq2, 'EdgeColor','none');
view(ax2,2)
xlim(ax2,[L2(1),L2(2)])
ylim(ax2,[L2(3),L2(4)])
clim(ax2, [gmin, gmax]);
axis(ax2,'tight')
xlabel(ax2,'$x$','Interpreter','latex')
ylabel(ax2,'$u_2,~y$','Interpreter','latex')
ax2.Position=[0.1 0.38 0.74 0.27];
set(ax2,'Color','none', ...
        'XTick',[], ...
        'Box','off');
ax2.FontSize=26;
ax2.YTick=[L2(3) L2(4)];
colormap(ax2,'jet')
ax2.YLabel.Units='normalized';
ax2.YLabel.Position=[-0.07 0.5 0];

ax1 = axes('Parent',fig);   % for layer-1 (bottom)
hold(ax1,'on');
u1_fig=surf(ax1, xq1, yq1, vq1, 'EdgeColor','none');
view(ax1,2)
xlim(ax1,[L1(1),L1(2)])
ylim(ax1,[L1(3),L1(4)])
clim(ax1, [gmin, gmax]);
axis(ax1,'tight')
xlabel(ax1,'$x$','Interpreter','latex')
ylabel(ax1,'$u_1,~y$','Interpreter','latex')
ax1.Position=[0.1 0.11 0.74 0.27];
ax1.FontSize=26;
ax1.YTick=[L1(3) L1(4)];
colormap(ax1,'jet')
ax1.YLabel.Units='normalized';
ax1.YLabel.Position=[-0.07 0.5 0];

set([ax1 ax2 ax3], 'CLim', [gmin,gmax])

colormap('jet'); 
cb_shared = colorbar('Position', [0.88, 0.11, 0.02, 0.81],'FontSize',24); 

saveID = 0;

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);

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

%%
if mod(ti,saveStride)==0 || ti==1

    saveID = saveID + 1;

    mobj0.tt(saveID,1) = t;

    mobj1.u11(:,saveID) = u1;
    mobj2.u22(:,saveID) = u2;
    mobj3.u33(:,saveID) = u3;

end

%% Plot and gif
    if mod(ti, drawperframe) == 1
        vq3 = griddata(p3(1,:),p3(2,:),u3,xq3,yq3);
        vq2 = griddata(p2(1,:),p2(2,:),u2,xq2,yq2);
        vq1 = griddata(p1(1,:),p1(2,:),u1,xq1,yq1);

        gmin = min([gmin; vq1(:); vq2(:); vq3(:)]);
        gmax = max([gmax; vq1(:); vq2(:); vq3(:)]);
        set([ax1 ax2 ax3], 'CLim', [gmin,gmax]);


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
            ax3.Title.String = sprintf('t = %.4f', t);
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
    if any(ti==jjj)
        fname = sprintf('%s%s_t_%3.2f.png', prefix, eval, t);
        saveas(fig, fname);
        fname = sprintf('%s%s_t_%3.2f.fig', prefix, eval, t);
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

%% kymograph values
    k_idx = k_idx + 1; % Increment the "packed" counter
    t_kymo(k_idx) = t;  % Save the current time
        
    % --- Kymograph: Averaging over X to see Y-profile changes ---
    % No transpose needed: mean(..., 2) gives (ny x 1)
    kx_u1_p(:, k_idx) = mean(vq1, 2, 'omitnan');
    kx_u2_p(:, k_idx) = mean(vq2, 2, 'omitnan');
    kx_u3_p(:, k_idx) = mean(vq3, 2, 'omitnan');

    % --- Kymograph: Averaging over Y to see X-profile changes ---
    % Transpose needed: mean(..., 1) gives (1 x nx), we need (nx x 1)
    ky_u1_p(:, k_idx) = mean(vq1, 1, 'omitnan').';
    ky_u2_p(:, k_idx) = mean(vq2, 1, 'omitnan').';
    ky_u3_p(:, k_idx) = mean(vq3, 1, 'omitnan').';
end

%% savemat
mobj0.p1 = p1;  mobj0.p2 = p2;  mobj0.p3 = p3;

mobj0.L1 = L1;  mobj0.L2 = L2;  mobj0.L3 = L3;

mobj0.xq1 = xq1;  mobj0.yq1 = yq1;
mobj0.xq2 = xq2;  mobj0.yq2 = yq2;
mobj0.xq3 = xq3;  mobj0.yq3 = yq3;

mobj0.saveStride = saveStride;
mobj0.dt = dt;
mobj0.T = T;
mobj0.nsave = saveID;

end
%% kymograph
save([prefix '_kymo_0.mat'], ...
     't_kymo','xq1','yq1','xq2','yq2', 'xq3', 'yq3', ...
     '-v7.3');
     
save([prefix '_kymo_kx_u1.mat'],'kx_u1_p','-v7.3')
save([prefix '_kymo_kx_u2.mat'],'kx_u2_p','-v7.3')
save([prefix '_kymo_kx_u3.mat'],'kx_u3_p','-v7.3')
save([prefix '_kymo_ky_u1.mat'],'ky_u1_p','-v7.3')
save([prefix '_kymo_ky_u2.mat'],'ky_u2_p','-v7.3')
save([prefix '_kymo_ky_u3.mat'],'ky_u3_p','-v7.3')

%
max_t=t;
max_y=max([ky_u1_p,ky_u2_p,ky_u3_p],[],'all');
min_y=min([ky_u1_p,ky_u2_p,ky_u3_p],[],'all');
skip=1;
% max_y=15;
% min_y=-15;
k_idxx=k_idx;

%
fig_ky1 = figure('Color','w');
imagesc(xq1(1,:), t_kymo(1:skip:k_idxx), ky_u1_p(:,1:skip:k_idxx).');
axis xy; 
ylabel('$t$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
xlim([min(xq1,[],'all'),max(xq1,[],'all')])
ylim([t_kymo(1),max_t])
clim([min_y,max_y])
set(gca,'FontSize',12)
ax1 = gca;
ax1.Position=[0.16,0.13,0.73,0.8];
colormap('jet'); 
cb1=colorbar;
cb1.Position=[0.92,0.13,0.02,0.8];

saveas(fig_ky1,[prefix, '_kymo_yavg_layer1.png']);
saveas(fig_ky1,[prefix, '_kymo_yavg_layer1.fig']);

%
fig_ky2 = figure('Color','w');
imagesc(xq2(1,:), t_kymo(1:skip:k_idxx), ky_u2_p(:,1:skip:k_idxx).');
axis xy; colormap('jet'); colorbar;
ylabel('$t$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
xlim([min(xq1,[],'all'),max(xq1,[],'all')])
ylim([t_kymo(1),max_t])
clim([min_y,max_y])
set(gca,'FontSize',12)
ax2 = gca;
ax2.Position=[0.16,0.13,0.73,0.8];
colormap('jet'); 
cb2=colorbar;
cb2.Position=[0.92,0.13,0.02,0.8];

saveas(fig_ky2,[prefix, '_kymo_yavg_layer2.png']);
saveas(fig_ky2,[prefix, '_kymo_yavg_layer2.fig']);

%
fig_ky3 = figure('Color','w');
imagesc(xq3(1,:), t_kymo(1:skip:k_idxx), ky_u3_p(:,1:skip:k_idxx).');
axis xy; colormap('jet'); colorbar;
ylabel('$t$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
xlim([min(xq1,[],'all'),max(xq1,[],'all')])
ylim([t_kymo(1),max_t])
clim([min_y,max_y])
set(gca,'FontSize',12)
ax3 = gca;
ax3.Position=[0.16,0.13,0.73,0.8];
colormap('jet'); 
cb3=colorbar;
cb3.Position=[0.92,0.13,0.02,0.8];

saveas(fig_ky3,[prefix, '_kymo_yavg_layer3.png']);
saveas(fig_ky3,[prefix, '_kymo_yavg_layer3.fig']);

%%
fprintf('\nDone!\n');
toc
diary OFF
