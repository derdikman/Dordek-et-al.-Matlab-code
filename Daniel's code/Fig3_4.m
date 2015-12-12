close all
clear all
clc

omega=1;% angular velocity
v=2.5; %velocity

comps=40; %number of PCA components to show
dims=[100 100]; %box dimensions (can be any number of dimensions)
D=length(dims);
sigma1=7.5; % positive Gaussian width
sigma2=2*sigma1; % negative Gaussian width
boundaries_set={0,'circular'}; % boundary condition options - look at imfilter boundary options for full options
boundaries=boundaries_set{2};  

win_size=sigma2*5; %Gaussian window width
gauss_diff_conv=@(z) imblur(z,sigma1*ones(1,D),win_size,D,boundaries) ...
            - imblur(z,sigma2*ones(1,D),win_size,D,boundaries) ;

T=1e6; 
agent_set={'white_noise','rand_walk'};
agent=agent_set{2};

batch=1e3;
if T<batch
    batch=T;
end

%     [C,S,lambda]=pca(Y,'NumComponents',comps);

reps=T/batch;
YY=0;
last_point=ceil(rand(1,D).*dims)';
last_angle=0;

for ii=1:reps
        
    switch agent
        case 'white_noise'
            X=randn([dims,batch]);
        case 'rand_walk'
            z=randn(1,batch);
            walk_D=last_angle+cumsum(omega*z,2);
            walk_D=mod(walk_D,2*pi);
            x_walk=mod(last_point(1)+round(v*cumsum(cos(walk_D))),dims(1))+1;
            y_walk=mod(last_point(2)+round(v*cumsum(sin(walk_D))),dims(2))+1;
            walk=[x_walk;y_walk];
            X=zeros([dims,batch]);
            linearInd = sub2ind(size(X), walk(1,:),walk(2,:),1:batch);
            X(linearInd )=1;
            last_point=walk(:,end);
            last_angle=walk_D(:,end);
    end
    Y=gauss_diff_conv(X);
    Y=reshape(Y,[prod(dims),batch])';
    YY=YY+Y'*Y; 
    disp(num2str(ii/reps))
end
    
[C,D]=eig(YY);
lambda=diag(D);
% end
lambda=lambda(end:-1:1);
C=fliplr(flipud(C));
last_ind=625;
C=C(:,1:last_ind);
lambda=lambda(1:last_ind);
%  save('pca_results','dims','sigma1','sigma2','boundaries','lambda','C','omega','v')
%%
figure(1)
lambda=lambda/max(lambda);
subplot(1,3,1);
plot(lambda,'.')
xlim([0,625])
xlabel('# of Eigenvalues')
ylabel('Eigenvalues')
subplot(1,3,2);
variance_explained=cumsum(lambda)/sum(lambda);
plot(variance_explained,'-','linewidth',3)
xlim([0,625])
xlabel('# of Eigenvalues')
ylabel('Variance Explained')
subplot(1,3,3);
plot(lambda(1:comps),'.')
xlabel('# of Eigenvalues')
ylabel('Eigenvalues')

% addpath(genpath('C:\Users\Daniel\Dropbox\Matlab\Auxiliary_functions'))
% Export2Folder('PCA_eigenvalues.tif',[])


figure(2)

for cc=1:comps
    subplot(4,comps/4,cc)
    im=reshape(C(:,cc),dims);
    imagesc(im);
end


% addpath(genpath('C:\Users\Daniel\Dropbox\Matlab\Auxiliary_functions'))
% Export2Folder('PCA_eigenvalues.tif',[])

figure(3)
comps=40;
lim=0.4;
for cc=1:comps
    subplot(4,comps/4,cc)
    im=reshape(C(:,cc),dims);
    fft_im=abs(fftshift(fft2(im)));
    imagesc(fft_im);
    xlim([lim 1-lim]*dims(1))
    ylim([lim 1-lim]*dims(2))
end

%%