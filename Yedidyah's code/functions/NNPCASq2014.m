function [ v] = NNPCASq2014( Conv,v)
%%%backup%%%%%%%%%%%%%
%%% Input_matrix = InputMatFixed; +  repmat(1e-3*randn(N1,1),1,L);
%%%%%%%%%%%%%%%%%%%%
%%% Input matrix needs to be of size pXn, n- samples, p - variables. The
%%% srcipt will transpose it, to get the data matrix X.

%%%%Initial conditions and %%%%%%%%%%%%%%%%%%%%%%%%%%
n = max(size(Conv)); %samples
p = min (size(Conv)); % variables
X = Conv;
if nargin <2
    v = rand(p,1) ;%v(t)
    v_prev = v;
else
    v(v<0) = 0;
end
T = 2e3;            %max time to run
fprintf('Starting calculation of Eigenvector\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:T
    %%%%%%THE CODE%%%%%%%%%%%%%%%%%%%%%%%%%%
    vp = v; vp(vp<0)=0;
    b = numel(vp(vp~=0)) / (sqrt(n) * norm(vp,2)+eps);
    v_prev1  = v;
    v=X*f(v,n)-b*f(v_prev,n);
    v_prev = v_prev1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%STOP conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (sum(isnan(v))>0)
        fprintf('Nan values!!!, breaking...\n')
        break
        end
    
    if(  norm(v_prev - v,'fro'))<1e-1
        fprintf('Algorithm converged, breaking...\n')
        break;
    end
    if (mod(t,10) ==0)
        fprintf('iteration # %d ',t)
        fprintf('Current converging distance is %d\n',norm(v_prev - v,'fro'));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(0.000001);  % for control purpose
end %end of loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Final cut offs of vectors%%%%%%%%%%%%%%%%%%
v(v<0) = 0; % cut off last negative parts of vector
v = v./(norm(v)+eps);
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
