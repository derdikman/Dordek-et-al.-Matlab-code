function [ NNPcaGridness,gridnessRegularPCA] = NNPCAGridness (InputMat,Res,ind_xS,ind_yS,N1,NN,v,placeCenters,a_std,b_std,c_std)
L=length(InputMat)-1; %time length simulated
%r_sd = (diff(r_s'))';
meanInputMat = mean(InputMat(:,1:L),2);
InputMatFixed = InputMat(:,1:L) - repmat(meanInputMat,1,L);
C = InputMatFixed*InputMatFixed';
Cov = (1/L)*C;
%% PCA
[EVec,signals1,EVal] = princomp(InputMat(:,1:L)');

minx = 0;
maxx = 10;
maxy = 10;
k=1;
beta = 0;
alpha1 = 1*10.^(1:1:5);
%find Gridness of NNPCA COV matrix
G= zeros(NN,NN,N1,k);

%%
      n_mat=mulMatCross(zeros(size(G,1)),zeros (size(G,1)));
for ii = 1:length(alpha1)
     M = hexWeights( N1,k,minx,maxx );
[ NNEVec, iterations ] = positiveEvec( Cov,alpha1(ii),beta,k,M);
        % figure; plot(NNEVec)


clear mat_cross
xx=linspace(0,maxx,NN);
yy=linspace(0,maxy,NN);
[R,C1] = ndgrid(xx, yy);

 for j=1:k
     
clear val; clear exponent;clear GC; clear exponent2;

    for i=1:length(placeCenters)
          exponent(:,:,i) = (a_std(i)*(R-placeCenters(1,i)).^2 +2*b_std(i)*(R-placeCenters(1,i)).*(C1-placeCenters(2,i)) + c_std(i)*(C1-placeCenters(2,i)).^2 );
                 val(:,:,i) = psiSat*(exp(-exponent(:,:,i)));
                if (v==1)
                    exponent(:,:,i) = (a_std(i)*(R-placeCenters(1,i)).^2 +2*b_std(i)*(R-placeCenters(1,i)).*(C1-placeCenters(2,i)) + c_std(i)*(C1-placeCenters(2,i)).^2 );
                    exponent2(:,:,i) = (a_std2(i)*(R-placeCenters(1,i)).^2 +2*b_std2(i)*(R-placeCenters(1,i)).*(C1-placeCenters(2,i)) + c_std2(i)*(C1-placeCenters(2,i)).^2 );
                    val(:,:,i) = psiSat*(exp(-exponent(:,:,i))) - psiSat/4*(exp(-exponent2(:,:,i)));
                G(:,:,j,ii) = G(:,:,j,ii) +  NNEVec(i,j)*val(:,:,i);
                else
                  G(:,:,j,ii) = G(:,:,j,ii) +  NNEVec(i,j)*val(:,:,i);
                end
    end

    mat_cross(:,:,j,ii)=Cross_Correlation(G(:,:,j,ii),G(:,:,j,ii),n_mat);

 end
 gridness(ii) = findGrindess(mat_cross(:,:,:,ii)) ;
 
end

%%
NNPcaGridness = max(findGrindess(mat_cross)) ;

%% Gridness of regular PCA
ind_xS(ind_xS==0)=1; %NEEDED???
ind_yS(ind_yS==0)=1;
totalActivity2 = zeros(Res,Res,N1);
occupancyTime2 = zeros(Res,Res);
for i=2:length(signals1)-1
    
    ind_x =ind_xS(i);
    ind_y =ind_yS(i);
    occupancyTime2(ind_x, ind_y) = occupancyTime2(ind_x, ind_y) + 1;
    
end

for j=1:size(signals1,2) %for all Evecs (or PCs)
    
    for i=2:length(signals1)-1
        
        ind_x =ind_xS(i); %where was the agent at the time (length of signals=T)
        ind_y =ind_yS(i);
        totalActivity2(ind_x, ind_y,j)  =totalActivity2(ind_x, ind_y,j)  + signals1(i,j);
    end
    
end

clear cross_map;
n_mat=mulMatCross(zeros(size(occupancyTime2)),zeros (size(occupancyTime2))) ;
for q=1:18 %up to 18 - the PCA outputs are still not degrgated.(check)
    K(:,:,q) =  (totalActivity2(:,:,q)./ (eps + occupancyTime2));
    cross_map(:,:,q)=Cross_Correlation(K(:,:,q),K(:,:,q),n_mat);
end


gridnessRegularPCA = findGrindess(cross_map)
fprintf('Finished calculating  Gridness. Moving to the next round...\n')

end
