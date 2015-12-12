clear all
close all
clc

% addpath(genpath('C:\Users\Daniel\Dropbox\Matlab\Auxiliary_functions'))
%% Plot params
SetDefaultGraphicSettings(1)
figure
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[2   2   12    12])
set(gcf, 'Position',[2   2   12  12])

title_font=12;
title_pos=[-0.1 1.05];
title_pos2=[-0.12 1.05];
%% FISTA params
sigma1=7.5; % positive Gaussian width
sigma2=sigma1*2; % negative Gaussian width
nCells=1; % number of modules
boundaries_set={0,'circular'}; % boundary condition options - look at imfilter boundary options for full options
boundaries=boundaries_set{2};  
iterations=2000; %number of Fista iterations (increase if we have convergence issues)
L=100; %size of box
%% 1D example
box_dims=L; %box dimensions (can be any number of dimensions)

[x1, Energy_array ] = RunFista(box_dims,[sigma1 sigma2],nCells,boundaries,iterations); %#ok

%% 2D example
box_dims=box_dims*[1 1]; %box dimensions (can be any number of dimensions)

[x2, Energy_array ] = RunFista(box_dims,[sigma1 sigma2],nCells,boundaries,iterations); %#ok

%% Check Analytical Predictions
L=box_dims(1);
func = @(x) -(exp(-0.5*sigma1^2*x.^2) - exp(-0.5*sigma2^2*x.^2));
k_star=fminbnd(func,0.1*2*pi/sigma2,10*2*pi/sigma1);
k_spacing=k_star/(2*pi/L);
grid_spacing=2*pi/k_star;

%%
second_color=[0.8500    0.3250    0.0980];
k_m=0.5;
subplot(3,1,1)
r=abs(fftshift(fft(x1)));
temp=0:(2*pi/L):pi;
k_lattice=[-temp(end:-1:1),temp(2:end-1)];
stem(k_lattice,r,'markersize',3)   
hold on 
line(k_star*[1,1],[0,10],'color',second_color)   
line(-k_star*[1,1],[0,10],'color',second_color)   
xlim([-k_m k_m])
title('A','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)

subplot(3,1,[2 3])
r=abs(fftshift(fft2(x2)));
imagesc(k_lattice,k_lattice,r)
% grid on
colorbar
hold on
x=-k_star:1e-5:k_star;
plot(x,sqrt(k_star^2-x.^2),'color',second_color)
hold on
plot(x,-sqrt(k_star^2-x.^2),'color',second_color)
hold on
for kx=k_lattice
    for ky=k_lattice
        hold all
        plot(kx,ky,'.w','markersize',1.5)   
    end
end
xlim([-k_m k_m])
ylim([-k_m k_m])
title('B','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold','fontsize',title_font)
% axis equal

% Export2Folder('Kspace.tif',[])
