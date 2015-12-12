function [ v] = NNPCA2014( Input_matrix,v)
%%%backup%%%%%%%%%%%%%
%%% Input_matrix = InputMatFixed; +  repmat(1e-3*randn(N1,1),1,L);
%%%%%%%%%%%%%%%%%%%%
%%% Input matrix needs to be of size pXn, n- samples, p - variables. The
%%% srcipt will transpose it, to get the data matrix X.

%%%%Initial conditions and %%%%%%%%%%%%%%%%%%%%%%%%%%
n = max(size(Input_matrix)); %samples
p = min (size(Input_matrix)); % variables
X = Input_matrix';
if nargin <2
    v = rand(p,1) ;%v(t)
   %  v = ones(p,1) ;%v(t)
    v(v<0)=0;
else
    v(v<0) = 0;
end
u= zeros(n,1) ;  %u(t)
T = 13e3;            %max time to run
fprintf('Starting calculation of Eigenvector\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:T
    %%%%%%THE CODE%%%%%%%%%%%%%%%%%%%%%%%%%%
    vp = v; vp(vp<0)=0;
    b = numel(vp(vp~=0)) / (sqrt(n) * norm(vp,2)+eps);
    u=X*f(v,n) - b*g(u,n);
    d = sqrt(n)/norm(u,2);
    v_prev  = v;
    v=X'*g(u,n) - d*f(v,n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STOP condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(  norm(v_prev - v,'fro'))<1e-1
        fprintf('Algorithm converged, breaking...\n')
        break;
    end
    if (t>2e3)
        fprintf('Iteration exceeded max iters allowed')
        break;
    end
    
    if (mod(t,50) ==0)
        fprintf('iteration # %d ',t)
        fprintf('Current converging distance is %d\n',norm(v_prev - v,'fro'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(0.000001);  % for control purpose
end %end of loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Final cut offs of vectors%%%%%%%%%%%%%%%%%%
v(v<0) = 0; % cut off last negative parts of vector
v = v./norm(v);
u(u<0) = 0;
u = u./norm(u);
fprintf('Done with calcaulations\n');
end %end of function

%% BACKUP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %crearting (v(t))+
%     vPlus = v;
%     vPlus(vPlus<0)=0;
%
%     b = numel(vPlus(vPlus~=0)) / (sqrt(n) * norm(vPlus,2));
%     %function of f or g, the only difference is the input type. added eps to prevent 0/0
%     f = @(x) sqrt(n)*x/(norm(x,2)+eps);
%     % calculating u
%     u = X*f(vPlus) - b*f(u_minus1);
%     %calculating v
%     d = sqrt(n)/norm(u,2); %use the new u
%    v_zero = X'*f(u) - d*f( vPlus);
%    u_minus1 = u;
%    v_save(:,t) = v_zero;
%%%%%%%%%%%%%%%%%%%%%%%
