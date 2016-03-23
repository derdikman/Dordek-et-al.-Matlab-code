% clear all
close all
clc

run_place_cell_parameters

%% Plot params
% load place_data_periodic_LargeBox


figure
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[2   2   12    12])
set(gcf, 'Position',[2   2   12  12])
orientation = squeeze(orientation);
scale = squeeze(scale);
gridness = squeeze(gridness);
title_pos=[-0.1 1.05];
title_font=12;

take = 1:13;
% take = 1:12;
sigma_vector = sigma_vector(take);
gridness = gridness(take,:);
orientation = orientation(take,:);
scale = scale(take,:);

[nscales niters] = size(gridness);

sigma_mat = repmat(sigma_vector',[1 niters]);

%%

scale_analytical=zeros(length(sigma_mat(:)),1);

for kk=1:length(sigma_mat(:))
    sigma1=sigma_mat(kk);
    sigma2=2*sigma1;

    func = @(x) -(exp(-0.5*sigma1^2*x.^2) - exp(-0.5*sigma2^2*x.^2));
%     func = @(x) -((2*normcdf(L/sigma1)-1)*exp(-0.5*sigma1^2*x.^2) - (2*normcdf(L/sigma2)-1)*exp(-0.5*sigma2^2*x.^2));
    k_star=fminbnd(func,0.1*2*pi/sigma2,10*2*pi/sigma1)
    scale_analytical(kk)=(2/sqrt(3))*2*pi/k_star;
end

reg_coeffs=[ones(length(sigma_mat(:)),1), sigma_mat(:)]\scale(:);
scale_fitted=sigma_mat(:)*reg_coeffs(2)+reg_coeffs(1);
figure;

%%
red_color=[0.8500    0.3250    0.0980];
subplot(3,2,1);
plot(sigma_mat(:),scale(:),'*',sigma_mat(:),scale_fitted(:),'-', 'linewidth',2);
hold all
plot(sigma_mat(:),scale_analytical(:),'-','color',red_color,'linewidth',2);
ylabel('scale');
lgnd=legend('simulation','linear fit','lower bound');
set(lgnd,'location','northwest')
title('A','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
xlim([sigma_vector(take(1)),sigma_vector(take(end))])

subplot(3,2,2);
plot(gridness(:),orientation(:),'*');
xlabel('gridness');
ylabel('orientation');
title('B','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
xlim([0.6,1.41])

subplot(3,2,3);
plot(sigma_mat(:),orientation(:),'*');
ylabel('orientation');
mean_orientation=mean(orientation,2);
hold all
plot(sigma_vector,mean_orientation,'-', 'linewidth',2);
title('C','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
xlim([sigma_vector(take(1)),sigma_vector(take(end))])
ylim([0 15])
grid on

subplot(3,2,4);
hist(orientation(:));
xlabel('orientation');
ylabel('count');
text(0.4,0.9,sprintf('median=%f',median(orientation(:))), 'Units', 'normalized');
title('D','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)


subplot(3,2,5);
plot(sigma_mat(:),gridness(:),'*');
xlabel('sigma');
ylabel('gridness');
title('E','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
xlim([sigma_vector(take(1)),sigma_vector(take(end))])


subplot(3,2,6)
hist(gridness(:));
xlabel('gridness');
ylabel('count');
text(0.1,0.9,sprintf('mean=%f',mean(gridness(:))), 'Units', 'normalized');
title('F','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold','fontsize',title_font)
xlim([0.6,1.41])

% addpath(genpath('C:\Users\Daniel\Dropbox\Matlab\Auxiliary_functions'))
% Export2Folder('SigmaScan_periodic.tif',[])
% Export2Folder('SigmaScan_0.tif',[])