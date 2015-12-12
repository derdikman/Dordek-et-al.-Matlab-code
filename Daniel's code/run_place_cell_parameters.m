function run_place_cell_parameters

N = 500; % size of box

sigma_vector = 2:0.5:10; % range of place cell widths to check
niters = 10; % run every scale a few times
name='place_data_periodic_LargeBox'
load_previous=0;

if load_previous
    load([name, '.mat'])
    start_i=ii;
else
    start_i=1;
end

for ii = start_i:length(sigma_vector)
    for j = 1:niters
        sigma1 = sigma_vector(ii)
%         N = sigma1*40;
        [gridness(:,ii,j),scale(:,ii,j),orientation(:,ii,j),rate_map{ii,j}] = grid_alignment(N,sigma1);
        save([name, '.mat'])
    end
end

disp('');
