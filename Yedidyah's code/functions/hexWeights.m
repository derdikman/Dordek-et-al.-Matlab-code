function [ J ] = hexWeights( N1,k,minx,maxx )

J = zeros(k,N1);
for g=1:k
    phi=15*rand;
    c=(1:3)';
    k=[cos(2*pi*c/3+phi), sin(2*pi*c/3+phi)];
    k_star=norm(k);
    k=k_star*k;
    
    x1=linspace(minx,maxx,sqrt(N1));
    psi_3=zeros(length(x1),length(x1));
    [ q, w] = ndgrid(x1,x1);
    noise =0;% .3* randn(length(x1)); %Gaussian noise
    for jj=1:3
        psi_3=cos( k(jj,1)*q+k(jj,2)*w+10*rand)+psi_3;
    end
    %Weights will be initialized as positive all!
      psi_3=(2/3)*psi_3+1+ noise;
       
%figure; imagesc(psi_3), colorbar
     J(g,:) = reshape(psi_3,1,N1);
     normFactor = repmat((sum(J(g,:,1))),1,N1);
     J(g,:) = 6*J(g,:,1)./normFactor;
%   %  
 end

J = J';
 J = J./repmat( sqrt(sum(J.^2,1)) , N1,1 ); 
end

