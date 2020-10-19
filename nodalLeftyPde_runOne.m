
%%
clearvars
close all

load('/Users/sapnachhabra/Desktop/CellTrackercd/model/NodalLeftySmad2/experimentalData/betaCatenin.mat');
bCat = experimental_bCat.radialAverageInterpolated;


%%%System Parameters

tic;
%% -------------------parameters
%[Nodal, Lefty, Lefty inhibitor, Wnt];
param.dC = [0.0008,0.0008,0.0000]; % diffusion constant
param.B1 = 0.1; param.K11_1 = 0.11; param.K11_2 = 1; param.n11 = 5;  param.kd1 = 0.0007; param.K21 = 0.4; %1;
param.K41_1 = 0.25; param.K41_2 = 1; param.n41 = 5;
param.B2 = 0.1;  param.K12_1 = 0.05; param.K12_2 = 1; param.n12 = 4; param.K32 = 0.01; param.kd2 = 0.001;
param.K13_1 = 0.06; param.K13_2 = 1.5;  param.n13 = 4;  param.kd3 = 0.0002;

%%
param.kdm = 10; % degradation of components outside colony
param.tmax = 1000;
param.dt = 1;

%%%initial conditions
param.ic(1) = 0.32;% 0.32 for initial phase activator
param.ic(2) = 2; %inhibitor, high enough to repress activator
param.ic(3) = 0.0001;

%% -------------------Simulation domain
colonyRadius = 3; %colony radius
edgeWidth = 0.5;

nx = 1024/4;
ny = 1024/4;

Mx = 3;
Lx = Mx*pi; %domain width
My = 3;
Ly = My*pi; %domain height
xgrid = linspace(-Lx, Lx, nx);
ygrid = linspace(-Ly, Ly, ny);
[X, Y] = meshgrid(xgrid,ygrid);
dx=2*Lx/nx; % effective discretization width x
dy=2*Ly/ny; %  effective discretization width y; for record keeping only

gridSpace = sqrt(X.^2 + Y.^2);
COL = 1-heaviside(gridSpace - colonyRadius);
[B,~] = bwboundaries(imbinarize(COL));

edge = COL - (1 - heaviside(gridSpace - (colonyRadius - edgeWidth)));
%%
%% ------------------initial conditions
u1_0 = ones(nx,ny).*(param.ic(1)).*COL; % nodal - uniformly distributed
%u1_0 = (param.ic(1)).*edge +  0.1.*ones(nx,ny).*(COL - edge);

u2_0 =  ones(nx,ny).*(param.ic(2)).*COL;
u3_0 =  rand(nx,ny).*(param.ic(3)).*COL; % lefty inhibitor - uniformly distributed

figure;
subplot(1,3,1); imagesc(u1_0);
subplot(1,3,2); imagesc(u2_0);
subplot(1,3,3); imagesc(u3_0);
%%
%% -------------------solve pde
% spectral solver part1
%%
kx = [[0:nx/2] [-nx/2+1: -1]]./Mx;
ky = [[0:ny/2] [-ny/2+1: -1]]./My;
%
nL=zeros(ny,nx);% negative Laplacian
for jj = 1:ny
    nL(jj,:) =  (ky(jj)^2+kx.^2);
end
%%differential operators
for ii = 1:numel((param.dC))
    diffOp(:,:,ii) = [1 + param.dt*(param.dC(ii)*nL)];
end
%%
%%

%%Define Nonlinearities
%% original equation
f1 = @(uu1,uu2,uu3,uu4) (param.K11_1).*(uu1.^param.n11)./((param.K11_2^param.n11 + uu1.^param.n11).*(1+param.K21.*uu2)) ...
    - (param.kd1).*uu1 + param.K41_1.*(uu4.^param.n41)./(param.K41_2.^param.n41+ uu4.^param.n41);
f2 = @(uu1,uu2,uu3,uu4) (param.K12_1).*(uu1.^param.n12)./(param.K12_2^param.n12 + uu1.^param.n12)...
    - (param.kd2).*uu2 - (param.K32).*uu3;
f3 = @(uu1,uu2,uu3,uu4) (param.K13_1).*(uu1.^param.n13)./(param.K13_2^param.n13 + uu1.^param.n13) - (param.kd3).*uu3 ;

%%
% spectral solver part2
%%Calculate fourier transforms of these initial data, v denote
%%fourier transformed variables
v1 = fft2(u1_0);
v2 = fft2(u2_0);
v3 = fft2(u3_0);

t = 0;


counter = 1;
%

%%%%Start the time stepper
while t< (param.tmax)
    t = t+(param.dt);
    u1 = max(0, real(ifft2(v1)));
    u2 = max(0, real(ifft2(v2)));
    u3 = max(0, real(ifft2(v3)));
    
    
    
    
    u4 = bCat(t*ones(nx,ny),gridSpace).*COL; % beta catenin
    u4(u4<0.6) = 0; 
    
    % --- uncomment for stopping wnt at 30h
    %         if t<30*(param.tmax)/45
    %             u4 = bCat(t*ones(nx,ny),gridSpace).*COL; % beta catenin
    %             u4(u4<0.6) = 0;
    %         else
    %             u4 = COL.*0;
    %         end
    
    
    
    
    if mod(round(t/(param.dt)),10)==0
        figure(4)
        
        subplot(2,2,1)
        imagesc(u1)
        title(['activator'])
        colorbar;
        %caxis([0 100])
        %hold on; plot(B{1}(:,1), B{1}(:,2), 'w-', 'lineWidth', 0.7);
        
        subplot(2,2,2)
        imagesc(u2)
        title(['inhibitor protein'])
        %caxis([0 1])
        colorbar
        %hold on; plot(B{1}(:,1), B{1}(:,2), 'w-', 'lineWidth', 0.7);
        
        subplot(2,2,3)
        imagesc(u3)
        title(['inhibitor of inhibitor'])
        colorbar
        hold on; plot(B{1}(:,1), B{1}(:,2), 'w-', 'lineWidth', 0.7);
        
        
        subplot(2,2,4);
        imagesc(u4); %imagesc(xgrid,xgrid,u4)
        title(['wnt'])
        colorbar; %caxis([0 1]);
        hold on; plot(B{1}(:,1), B{1}(:,2), 'w-', 'lineWidth', 0.7);
        
        
        
        figure(5); %plot a cross-section of u and v concentration
        subplot(3,1,1)
        plot(xgrid,real(u1(end/2,:))); xlim([-colonyRadius-0.5 +colonyRadius+0.5]); %ylim([0 1.8]);
        subplot(3,1,2)
        plot(xgrid,real(u2(end/2,:))); xlim([-colonyRadius-0.5 +colonyRadius+0.5]); %ylim([0 1]);
        subplot(3,1,3)
        plot(xgrid,real(u3(end/2,:))); xlim([-colonyRadius-0.5 +colonyRadius+0.5]);
        title([' t = ' num2str(t) ]);
    end
    
    
    v1 =( v1 + (param.dt).*fft2(COL.*f1(u1,u2,u3,u4) - (1-COL).*u1.*(param.kdm)   ) )./diffOp(:,:,1);
    v2 =( v2 + (param.dt).*fft2(COL.*f2(u1,u2,u3,u4) - (1-COL).*u2.*(param.kdm)  ) )./diffOp(:,:,2);
    v3 =( v3 + (param.dt).*fft2(COL.*f3(u1,u2,u3,u4) - (1-COL).*u3.*(param.kdm)  ) )./diffOp(:,:,3);
    
    
end
%%
toc;