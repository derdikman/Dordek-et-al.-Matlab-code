function [gridness,scale,orientation,rate_maps] = grid_alignment(N,sigma1)
%
% Use Fista algorithm to calc Positive-PCA and then ged the alignment of
% the grid to the closest wall
tic

if ~exist('N') || isempty(N)
    %N = 200;
    N=100;
end

box_dims=[N N]; %box size

if ~exist('sigma1') || isempty(sigma1)
    sigma1=10; % positive Gaussian width
    %sigma1 = 5;
end

sigma2=2*sigma1; % negative Gaussian width
boundaries_set={0,'circular','replicate'}; % boundary condition options - look at imfilter boundary options for full options
boundaries=boundaries_set{1};

nCells = 1; 
iterations = 1000;
mask= 1; %rand([box_dims,nModules])<0.5;% either 1 or a boolean matrix of size [dims, nModules] - indicating locations of place cells

% figure(666);
[rate_maps, Energy_array ] = RunFista(box_dims,[sigma1 sigma2],nCells,boundaries,iterations);
%[rate_maps, Energy_array ] = RunFista(box_dims,[sigma1 sigma2],nModules,lambda,boundaries,mask,iterations);

% title(sigma1);
% drawnow;
% 
% figure(777);
% hold off;
for i = 1:nCells
    [gridness(i),scale(i),orientation(i)] = calc_grid_alignment(rate_maps(:,:,i));
%     drawnow;
%     pause(1);
end

save grid_alignment_data
toc

% figure(888);
subplot(2,1,1);
plot(scale,'.');
subplot(2,1,2);
plot(gridness,'.');
disp('');

% figure(999);
% plot(scale(gridness > 0.6),'.');

disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gridness,scale,orientation] = calc_grid_alignment(rate_map)

% calc the spatial autocorrelation plot of the rate map

autocorr_mat = Cross_Correlation(rate_map);

%fft_mat=abs(fftshift(fft2(rate_map)));

% imagesc(autocorr_mat);
% axis equal; axis off;
% hold on;

gridness = findGridness(autocorr_mat);

hsize = [9 9]; sigma = 3;
h = fspecial('gaussian', hsize, sigma);
autocorr_mat = imfilter(autocorr_mat,h);

peaks = find_peaks(autocorr_mat);
% plot(peaks(:,2),peaks(:,1),'.');

[orientation,scale] = find_orientation_from_peaks(peaks,autocorr_mat);

% title(sprintf('gridness = %f  scale = %f orientation = %f', ...
%               gridness, scale, orientation));

disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peaks = find_peaks(mat)

mat_center = mat(2:end-1,2:end-1);

mat_up = mat(1:end-2,2:end-1);
mat_down = mat(3:end,2:end-1);
mat_left = mat(2:end-1,1:end-2);
mat_right = mat(2:end-1,3:end);

peak_mat = mat_center > mat_up & mat_center > mat_down & ...
           mat_center > mat_left & mat_center > mat_right;
       
[i,j] = find(peak_mat);

peaks = [i+1 j+1];     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [orientation,scale] = find_orientation_from_peaks(peaks,autocorr_mat);
%
% find orientation of grid (relative to two walls) and grid scale,
% from six closest peaks to center in
% autocorrelation plot
%
[N M] = size(autocorr_mat);

% find central-peak in peaks (which is peak closest to the one you think
% should be in the center)
%
central_peak = [(N+1)/2 (M+1)/2];

dist_to_central_peak = sqrt( (peaks(:,1)-central_peak(1)).^2 + ...
                             (peaks(:,2)-central_peak(2)).^2 );
[min_dist, min_ind] = min(dist_to_central_peak);
    
central_peak = peaks(min_ind,:);
dist_to_central_peak = sqrt( (peaks(:,1)-central_peak(1)).^2 + ...
                             (peaks(:,2)-central_peak(2)).^2 );

%
% find six closest peaks to central peak
%

[sorted_dists,dist_inds] = sort(dist_to_central_peak);

closest_peaks = peaks(dist_inds(2:7),:);

%
% get grid scale
%

scale = median(sorted_dists(2:7));

% 
% get grid orientations for 6 closest peaks
%

peaks_vector(:,1) = closest_peaks(:,1) - central_peak(1);
peaks_vector(:,2) = closest_peaks(:,2) - central_peak(2);

% 
% get cosine of angles to vectors [0 1] and [1 0], representing walls
%

angles = [];
for i = 1:size(peaks_vector,1);
    angles(end+1) = subspace(peaks_vector(i,:)',[0 1]');
end
for i = 1:size(peaks_vector,1);
    angles(end+1) = subspace(peaks_vector(i,:)',[1 0]');
end

%disp(rad2deg(angles));

orientation = min(rad2deg(angles))

disp('');
