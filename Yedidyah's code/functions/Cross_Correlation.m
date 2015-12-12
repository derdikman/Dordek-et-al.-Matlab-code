function out_mat=Cross_Correlation(mat1,mat2,n)
% Input:mat1, mat2 assumed to be the same size. n: matrix holding the number
% of cells joining multiplication in both matrices per time unit.
%output: correlation matrix
% if nargin<2
%     mat2=mat1;
%     n=mulMatCross(zeros(length(mat1)),zeros(length(mat2)));
%   elseif nargin<3
%     n=mulMatCross(zeros(length(mat1)),zeros(length(mat2)));
% end
%%%%%%%%%%%
%get rid of ?NaNs
mat1(isnan(mat1)) = 0;
mat2(isnan(mat2)) = 0;

%%%%%%%%%%%%

Zmat=ones(2*size(mat2)-1); %used for a shift right & down - zeroing relavent culmn& row
Zmat(1,:)=0;
Zmat(:,1)=0;
onesMat=ones(size(mat2)); %for interior summation
%basically doing this: cov(x,y)/(sigma(x)*sigma(y). see http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
sumXY=circshift(xcorr2(mat1,mat2),[1 1]).*Zmat;
sumX=circshift(xcorr2(mat1,onesMat),[1 1]).*Zmat;
sumY=circshift(xcorr2(onesMat,mat2),[1 1]).*Zmat;
xSq=circshift(xcorr2(mat1.^2,onesMat),[1 1]).*Zmat;
ySq=circshift(xcorr2(onesMat,mat2.^2),[1 1]).*Zmat;
out_mat=(n.*sumXY-sumX.*sumY)./( (sqrt( n.*xSq- (sumX).^2  )  .* sqrt((n+1).*ySq-(sumY).^2)) +eps);

end