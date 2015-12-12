%%
clear all
close all
%%
k = N1;
run_len =5e5;
%%
X= randn(500,k);
S=size(X,1);
meanX = mean(X,1);
 XFixed =X - repmat(meanX,S,1);
C1 = XFixed'*XFixed;
%%
C1 = Cov(:,p);
%w = zeros(run_len,k);
% w2=zeros(1,k);
% w2 = rand(1,k);

w2 = 1/k.*ones(1,k);  %choose random values for init

iter =2;
alpha = 1e-2; %advance size
while ( (iter<run_len) ) %&& (sum(alpha* (2*C*w2')') > 1e-4)  )
    %advance rule
    w2 =   w2 + alpha* (2*C1*w2')';
    %take max(0,Wk) to keep nonnegativity
    w2(w2<0)=0;
      %make sure norm=1
      w2 =   w2./norm(w2);
      %for later analyze
      %sum_w(iter)=sum(alpha* (2*C*w2')'); %sum of diffs
      %advance counter
      iter = iter+1
     %  W(iter,:)=w2;
end

