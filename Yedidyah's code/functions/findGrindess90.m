function [gridness, gridnessSTD]=findGrindess90(cross_map)
%usage
auto_cross_original=cross_map;
rot_auto_corr=zeros(size(cross_map,1),size(cross_map,2),size(cross_map,3),3);
for m=1:3
    rot_auto_corr(:,:,:,m)=imrotate(auto_cross_original,m*45,'bicubic','crop');
end

Xc = size(cross_map,1)/2;
Yc=size(cross_map,2)/2;

%find the inner circle of the autocorrelation's radiuss'
% meanMiddleRadi=0
% y2 = auto_cross_original(:,:,1);
%  r1=5;
%  r2=15;
%  [~, radii,~] = imfindcircles(y2,[r1 r2]);
%
%   for e=1:size(rot_auto_corr,3)
%  y2 = auto_cross_original(:,:,e);
%   r1=5;
%   r2=15;
%   [~, radii,~] = imfindcircles(y2,[r1 r2]);
%  if (isempty(radii))
%       meanMiddleRadi(e) = mean(meanMiddleRadi);
%   else
%       meanMiddleRadi(e)=mean(radii);
%   end
%
%   end


%start checking from the end.
radi_end = size(auto_cross_original,1)/1.1;
%take the mean of the radii found as the middle one. (I know. but do u have
%a better sugestion?)% radi_start = mean(radii);
radi_start = 6;

[Y,X] = ndgrid(1:1:Xc*2, 1:1:Yc*2);
radi=radi_end:-1:radi_start+40;

gridness=zeros(size(rot_auto_corr,3),size(rot_auto_corr,4));
best_ring=zeros(size(rot_auto_corr,3),size(rot_auto_corr,4));

for i=1:size(rot_auto_corr,3) %per autocorr matrix (grid cell unit correlation)
    
    tmp_gridness = zeros(length(radi),size(rot_auto_corr,4));
    for j=1:length(radi) %per radii
        I=zeros(size(Y));
        I( ( ((X-Xc).^2 + (Y-Yc).^2) < radi(j)^2 ) & (((X-Xc).^2 + (Y-Yc).^2) > radi_start^2) ) =1;
        ring_orignal=I.*auto_cross_original(:,:,i);
        
        %         I(((X-Xc).^2 + (Y-Yc).^2) > radi_start^2)=1;
           %      figure; imagesc(ring_orignal)
        
        
        %    ring_orignal(ring_orignal==0)=NaN;
        
        for m=1:size(rot_auto_corr,4)   %per rotation (we have 3)
            
            ring_cross=I.*rot_auto_corr(:,:,i,m);
            %ring_cross(ring_cross==0)=NaN;
            % tmp_gridness(j,m)=PointCorr(ring_cross,ring_orignal);
            tmp_gridness(j,m)= [0 1]*(corrcoef(ring_cross,ring_orignal))*[1 0]';
            %             if (mod(m,2)==0)   %for even rotations (60,120) the best gridness is the biggest
            %                 [gridness(i,m) ,best_idx]= max(tmp_gridness(:,m));
            %                 best_ring(i,m) = radi(best_idx);
            %             else                %for odd rotations (60,120) the best gridness is the smallest (looking for anti-corr)
            %                 [gridness(i,m) ,best_idx ]= min(tmp_gridness(:,m));
            %                 best_ring(i,m) = radi(best_idx);
            %             end
        end
    end %now check what is the best radi and gridness per rotation
    
    [gridness(i,[2 ]) ,best_idx]= max(tmp_gridness(:,[2]));
    best_ring(i,[2 ]) = radi(best_idx);
    [gridness(i,[1 3 ]) ,best_idx]= min(tmp_gridness(:,[1 3]));
    best_ring(i,[1 3 ] )= radi(best_idx);
    
end


final_gridness2=1*(gridness(:,2))-1/2*(gridness(:,1)+gridness(:,3));

M = mean(final_gridness2);
S = std(final_gridness2);
% mean1=mean(final_gridness1);
% mean2=mean(final_gridness2);
    [gridness, idx] = min(M);
    gridnessSTD =S(idx);

end
