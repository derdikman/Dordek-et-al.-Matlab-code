function [Struct] = PCAEvecCalc(Struct,r_saved,NumberOfPC)
%This function returns the regular and nonnegative PC of the data sent to
%it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs: r_saved - data [features,samples], NumberOfPC - #of PC needed to caclulate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Calculating PCA...\n')

% get the size of the data matrix. make sure the dimensions are alighned. 
[M,N] = size(r_saved);
if N<M
    r_saved = transpose(r_saved);
end
% Calculate the regular unconstrained PC of the data using eig function
 [~,PCAEvec,Evals] = pca1(r_saved);
 Struct.PCAEvec=PCAEvec;
 % save for output only # of PC required by used
 %PCAEvec = PCAEvec(:,1:NumberOfPC);
 L=length(r_saved');
 %%%%%%%%%%
%print images from paper
     C = r_saved*r_saved';
    Cov = (1/L)*C; 
    fprintf('printing images from paper...\n')
      figure;
    imagesc(r_saved(:,1:1e4)), colormap jet
    title('figure 1C')
    figure;
    imagesc(Cov), colormap jet
  title('figure 1D')
figure;
plot(Evals,'b.'), colormap jet
    title('figure 4A')
    figure;
    plot(cumsum(Evals)./sum(Evals),'b.'), colormap jet
        title('figure 4B')
  
%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Non-negative PCA
%  fprintf('Calculating Non-Negative PCA...\n')
% 
% for z=1:NumberOfPC %per PC 
%  % subtract the mean from every feature
% meanInputMat = mean(r_saved(:,1:L),2);
% InputMatFixed = r_saved(:,1:L) - repmat(meanInputMat,1,L);
% 
% 
% 
% 
% % send the fixed data to the Monatrani algorithm, see (Montanari A, Richard E. Non-negative principal component analysis: 
% %Message passing algorithms and sharp asymptotics. arXiv preprint arXiv:14064775. 2014.
% 
% PCANNEvec(:,z)=NNPCA2014( InputMatFixed);
% % check if there's Nan's in the eigenvectors!
% if sum( isnan(PCANNEvec(:,z))) >1
%     fprintf('NaNs!! breaking...\n');
%     break;
% end
% % caclulate the corresponding eigenvalue
% alpha1(:,z) = PCANNEvec(:,z)'*(InputMatFixed*InputMatFixed' )*PCANNEvec(:,z);
% % subtract the projection of the data on the 1st eigenvector from the data,
% % and reuse the "new data" for the next PC calc.
% InputMatFixed = InputMatFixed - PCANNEvec(:,z) * (PCANNEvec(:,z)'*InputMatFixed);
% 
% end
%Struct.PCANNEvec=PCANNEvec;
end