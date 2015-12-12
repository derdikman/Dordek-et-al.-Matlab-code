function [unitMaps] = Sharpen( totalActivity,occupancyTime)


OccupancyMat=repmat(occupancyTime,[1,1,size(totalActivity,3)]);
K = (totalActivity./ (eps + OccupancyMat));
K(K<0)=0;
%LPF
 X_aux = 1:1:size(totalActivity,1)+10;
 [XX, YY]=ndgrid(X_aux,X_aux);
[B,A] = butter(10,0.8,'low') ;
w=5;
h = 1/w^2*ones(w,w);
clear Z1;
for i=1:size(K,3)
    A1 =filter2(B,A,K(:,:,i));
    A2=imfilter(A1,h);
    Z1(:,:,i) =  A2./max(max( A2));

% q=27;
% figure; imagesc(K(:,:,i)),colorbar
% figure; imagesc(A1),colorbar
% figure; imagesc(A2),colorbar
%  figure; imagesc(Z1(:,:,i)),colorbar
%  

    R=20;
        thresh = 10;

    p= FastPeakFind( Z1(:,:,i),thresh,ones(R)./R^2);
%    imagesc(Z1(:,:,i)); hold on
%    plot(p(1:2:end),p(2:2:end),'b+'), colorbar
   clear stack
   for j=1:2:length(p)
       if (Z1(p(j+1),p(j),i)>0.3)
           stack(j) = 1;
           stack(j+1) = 1;
       else
           stack(j) = 0;
           stack(j+1) = 0;
       end
   end
       
   p2 =  p.*stack';
   p3=p2(p2>0);
   
   centers =[ p3(1:2:end)';p3(2:2:end)']';
radii =R.*ones(size(centers,1),1) ;

%

%
 I=zeros(size(YY));
 for f=1:size(centers,1)
     I( ( ( (XX-centers(f,1)).^2 + (YY-centers(f,2)).^2) < radii(f)^2 ))=1;
 end
%  figure; imagesc(I')
 %y2 = Z1(:,:,i);
%
 unitMaps(:,:,i)=(I'.*Z1(:,:,i));
 
end
%  figure; imagesc(unitMaps(:,:,i))
   
 %
%y1 = K(:,:,1);
 %M=Cross_Correlation(unitMaps(:,:,i),unitMaps(:,:,i),n_mat1);
%  figure; imagesc(M), colorbar
%   [gridness, gridnessSTD]=findGrindess(M)
  
%unitMaps  = y3;


% 2 times averging
% w=25; %chose a window size of 15. can be changed...
% h = 1/w^2*ones(w,w);
% Y=imfilter(Z1(:,:,q),h);
% figure; imagesc(Y), colorbar
% unitMaps =  imfilter(Y,h);
%recalculating n_mat for later usage.. if needed
%n_mat1=mulMatCross(zeros(size(totalActivity,1)+10),zeros (size(totalActivity,1)+10)) ;

%% Debug & future
% L= length(Z1);
% NFFT = 2^nextpow2(L);
% Fs=4e3;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% Z2 =abs (fftshift(fft2(Z1)));
% figure;
% imagesc(Z2)
% 
% figure;
% plot(f,2*abs(Z2(1:NFFT/2+1))) 


% K2=unitMaps(:,:,1);
% figure; imagesc(K2)
% M=Cross_Correlation(K2,K2,n_mat);
% figure; imagesc(M), colorbar
% figure; imagesc(Z1)


% y = filter2(h,Z1);
%  figure; imagesc(y), colorbar
%  title(['w ',num2str(w)]);



% y2 = filter2(h,y);
% figure; imagesc(y2), colorbar

% M=Cross_Correlation(y2,y2,n_mat);
% figure; imagesc(M), colorbar
% [gridness, gridnessSTD]=findGrindess(M);
%
% r1=10;
% r2=25;
% [centers, radii, metric] = imfindcircles(y2,[r1 r2]);
% % centersStrong5 = centers(1:5,:);
% % radiiStrong5 = radii(1:5);
% % metricStrong5 = metric(1:5);
% viscircles(centers, radii,'EdgeColor','b');
% % [centers(4,1),centers(5,1)] = deal(centers(5,1),centers(4,1));
% % [centers(4,2),centers(5,2)] = deal(centers(5,2),centers(4,2));

   
% D = pdist(centers);
% Z = squareform(D)+11e3*eye((size(centers,1)),(size(centers,1)));
% [~, idx] = find(Z(:)<20);%recheck. also maybe if ther's 2 circles - take the one of largest radiuos.
% %check a beeter way to find intersecting circles
% [ind_i,~] = ind2sub(size(Z),idx);
%
%
% centersFixed = centers;
% radiiFixed = radii;
% centersFixed(ind_i,:)=[];
% radiiFixed(ind_i,:)=[];
% viscircles(centersFixed, radiiFixed,'EdgeColor','b');
%
%
%  centers=centersFixed ;
% radii=radiiFixed;
% centers =[ p(1:2:end)';p(2:2:end)']';
% radii =20.*ones(size(centers,1),1) ;
% 
% %
%  X_aux = 1:1:141;
%  [XX, YY]=ndgrid(X_aux,X_aux);
% %
%  I=zeros(size(YY));
%  for i=1:length(centers), I( ( ((XX-centers(i,1)).^2 + (YY-centers(i,2)).^2) < radii(i)^2 ))=1;end
%  figure; imagesc(I')
%  y2 = Z1(:,:,q);
% %
%  y3=(I'.*y2);
%  figure; imagesc(y3)
% % r1=10;
% % r2=30;
% % [centers, radii, metric] = imfindcircles(y3,[r1 r2]);
% % viscircles(centers, radii,'EdgeColor','b');
% %
% y1 = K(:,:,1);
%  M=Cross_Correlation(y3,y3,n_mat1);
%  figure; imagesc(M), colorbar

end
