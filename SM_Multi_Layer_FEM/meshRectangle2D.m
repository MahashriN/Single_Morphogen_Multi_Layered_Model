clear all; close all; clc;
save_file=1; % 0-only see meshes; 1-save the mesh files

folder='D:\MATLAB Output\SM_Multi_Layer_FEM\meshes\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];

H=0.1;    
L=20;
dy=H/10; 
dx=L/250;
N=2;
pf=['N=',num2str(N),'\','L=',num2str(L),'\','H=',num2str(H),'\'];

%%
if N>2
L3=[0,L,2*H,3*H];
[p3,t3,p3_1DB,t3_1DB]=Rectangle2D(L3,dx,dy,50,'edge',[1,3],[prefix,'_3']);
bottom_p3 = p3_1DB(1,:);
bottom_t3 = t3_1DB(1:2,:);
top_p3    = p3_1DB(2,:);
top_t3    = t3_1DB(3:4,:);
if save_file
save([pf,'p3_2D.mat'],'-mat','p3')
save([pf,'t3_2D.mat'],'-mat','t3')
save([pf,'bottom_p3_1D.mat'],'-mat','bottom_p3')
save([pf,'bottom_t3_1D.mat'],'-mat','bottom_t3')
save([pf,'top_p3_1D.mat'],'-mat','top_p3')
save([pf,'top_t3_1D.mat'],'-mat','top_t3')
end
end
%%
if N>1
L2=[0,L,H,2*H];
[p2,t2,p2_1DB,t2_1DB]=Rectangle2D(L2,dx,dy,50,'edge',[1,3],[prefix,'_2']);
bottom_p2 = p2_1DB(1,:);
bottom_t2 = t2_1DB(1:2,:);
top_p2    = p2_1DB(2,:);
top_t2    = t2_1DB(3:4,:);
if save_file
save([pf,'p2_2D.mat'],'-mat','p2')
save([pf,'t2_2D.mat'],'-mat','t2')
save([pf,'bottom_p2_1D.mat'],'-mat','bottom_p2')
save([pf,'bottom_t2_1D.mat'],'-mat','bottom_t2')
save([pf,'top_p2_1D.mat'],'-mat','top_p2')
save([pf,'top_t2_1D.mat'],'-mat','top_t2')
end
end
%%
if N>0
L1=[0,L,0,H];
[p1,t1,p1_1DB,t1_1DB]=Rectangle2D(L1,dx,dy,50,'edge',[1,3],[prefix,'_1']);
bottom_p1 = p1_1DB(1,:);
bottom_t1 = t1_1DB(1:2,:);
top_p1    = p1_1DB(2,:);
top_t1    = t1_1DB(3:4,:);
if save_file
save([pf,'p1_2D.mat'],'-mat','p1')
save([pf,'t1_2D.mat'],'-mat','t1')
save([pf,'bottom_p1_1D.mat'],'-mat','bottom_p1')
save([pf,'bottom_t1_1D.mat'],'-mat','bottom_t1')
save([pf,'top_p1_1D.mat'],'-mat','top_p1')
save([pf,'top_t1_1D.mat'],'-mat','top_t1')
end
end
%%
function [p,t,bp,bt]=Rectangle2D(L,dx,dy,d,RegionType,RegionID,prefix)
Lxa=L(1); Lxb=L(2);
Lya=L(3); Lyb=L(4);
model=createpde();
scale = dx/dy;
R = [3,4,Lxa,Lxb,Lxb,Lxa, Lya*scale,Lya*scale,Lyb*scale,Lyb*scale]';
g = decsg(R);
geom = geometryFromEdges(model, g);
mesh = generateMesh(model,'Hmax',dx,"GeometricOrder","linear");
p=mesh.Nodes;
p(2,:) = p(2,:)/scale;
t=mesh.Elements; 
for i=1:length(RegionID)
bp(i,:) = findNodes(mesh,"region",RegionType,RegionID(i));
size(bp);
if RegionType=="Edge" || RegionType=="edge"
    if RegionID(i)==3
        bp(i,:)=flip(bp(i,:));
    elseif RegionID(i)==1
        bp(i,:)=bp(i,:);
    else
        error("Boundary edge is not top or bottom")
    end
else
    error("RegionType is not a string or not 'Edge' or 'edge'")
end
bt(2*i-1:2*i,:)=[bp(i,1:end-1);bp(i,2:end)];
end
fig=figure;
subplot(1,3,1)
pdegplot(model,EdgeLabels="on")
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
subplot(1,3,2)
pdemesh(model,NodeLabels="on",ElementLabels="on")
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
subplot(1,3,3)
hold on;
for k = 1:size(bt,1)
    idx = bt(k,:);
    plot(p(1,idx), p(2,idx), 'g*-', 'LineWidth', 1);
end
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
saveas(fig,[prefix,'_mesh.png']);
saveas(fig,[prefix,'_mesh.fig']);
end
