clear all
close all
clc


%% Plot params
SetDefaultGraphicSettings(1)
figure
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[2   2   19    12])
set(gcf, 'Position',[2   2   19  12])
title_font=12;
title_pos=[-0.05 1.05];
title_pos2=[-0.025 1.025];
%% 
sigma1=7.5; % positive Gaussian width
sigma2=2*sigma1; % negative Gaussian width
L=100;
c=5;

p_k = @(x) -(exp(-0.5*sigma1^2*x.^2) - exp(-0.5*sigma2^2*x.^2));
p_x = @(x) (exp(-0.5*x.^2/sigma1^2)/sigma1 - exp(-0.5*x.^2/sigma2^2)/sigma2);

x_lim=sigma2*3;
x=(-1:0.01:1)*x_lim;
k_lim=3*pi/sigma2;
k=k_lim*(-1:0.01:1);
k_lattice=0:(2*pi/L):k_lim;
k_lattice=[-k_lattice(end:-1:2), k_lattice];
k_star=fminbnd(p_k ,0.04,4);


subplot(2,3,1)
plot(x,p_x(x))
xlabel('$x$','Interpreter','Latex')
ylabel('$r\left(x\right)$','Interpreter','Latex')
xlim([-x_lim,x_lim])
title('A','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
subplot(2,3,4)
plot(k,-p_k(k))
hold all
plot(k_star*[-1,1],-p_k(k_star*[-1,1]),'x')
hold all
plot(k_lattice,-p_k(k_lattice),'ok','markersize',5)
xlim([-k_lim,k_lim])
ylim([0,0.5])
xlabel('$k$','Interpreter','latex')
ylabel('$\hat{r}\, \left(k\right)$','Interpreter','latex')
title('B','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)

subplot(2,3,[2 3 5 6])
x=-k_star:1e-4:k_star;
plot(x,sqrt(k_star^2-x.^2),'color',[0.8500    0.3250    0.0980])
hold on
plot(x,-sqrt(k_star^2-x.^2),'color',[0.8500    0.3250    0.0980])
k_lim=k_star*1.5;
k_list=0:(2*pi/L):k_lim;
k_list=[-k_list(end:-1:2), k_list];
for kx=k_list
    for ky=k_list
        hold all
        plot(kx,ky,'ok','markersize',5)
    end
end

k_space=0:0.001:k_lim;
k_space=[-k_space(end:-1:2), k_space];
[KX,KY]=meshgrid(k_space,k_space);
PK=-p_k(sqrt(KX.^2+KY.^2));
hold all
contour(KX,KY,PK,10)
% colormap('jet')
colorbar


title('C','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold','fontsize',title_font)
% xlabel('$k_x$','Interpreter','latex')
% ylabel('$k_y$','Interpreter','latex')

%% Degeneracty

k_groups_cell=cell(0,1);
k_group=zeros(4,2);
k_group(:,1)=[2,-2,0,0]';
k_group(:,2)=[0,0,2,-2]';
k_groups_cell{end+1}=k_group*2*pi/L;

k_group=zeros(8,2);
k_group(:,1)=[2,1,-2,-1,-2,-1,2,1]';
k_group(:,2)=[1,2,1,2,-1,-2,-1,-2]';
k_groups_cell{end+1}=k_group*2*pi/L;

k_group=zeros(4,2);
k_group(:,1)=[-2,2,-2,2]';
k_group(:,2)=[2,2,-2,-2]';
k_groups_cell{end+1}=k_group*2*pi/L;

k_group=zeros(4,2);
k_group(:,1)=[1,-1,1,-1]';
k_group(:,2)=[1,1,-1,-1]';
k_groups_cell{end+1}=k_group*2*pi/L;


K=length(k_groups_cell);
goodness=zeros(K,1);
for kk=1:K
    k_group=k_groups_cell{kk};
    goodness(kk)=sum(abs(k_group(:,1).^2+k_group(:,2).^2-(k_star)^2))
    for ii=1:size(k_group,1)
        hold all
        letter=[char((kk-1)+'A'),num2str(ii)] ;
        text(k_group(ii,1)+0.005,k_group(ii,2),letter,'color','b','fontweight','bold')
    end
end

% print('LatticePoints',gcf,'-dpdf')

% set(gca,'xtick',[],'ytick',[])
%%

% Export2Folder('TuningCurve.tif',[])