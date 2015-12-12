function [ u, t] = positiveEvec( Cov,alpha1,beta,k,EVec)
% close all
% clear variables
% clc

%%
%init w0
% CStruct = load('Cov.mat');
% Cov = CStruct.Cov;
%  M = 196;
%  C = Cov(1:M,1:M) ;%+abs(min(min(Cov(1:M,1:M)))); %dosnt work on negative??
% %Cov = rand(100);
% % p = randperm(size(C,1)); %create randpermutation of the columns of Cov
% % C =Cov(p,:);

C = Cov;
% % figure; imagesc(C), colorbar
%alpha1 = 1e7; %try for example for 1e1 and 1e5 and see the difference!
d = size(C,1) ; %size of each vector
%k=5; %number of rows (vectors)
u = EVec(:,1:k);%ADD NOISE?
u(u<0)=0;
u=u./repmat( sqrt(sum(u.^2,1)) , d,1 ); 
%u=rand(d,k);          %initial guess for u!!! how does it affect it?
%u=u./repmat( sqrt(sum(u.^2,1)) , d,1 ); %normalize the columns of 'u' to unit norm!!

u_prev=zeros(size(u));
%beta = 0; %for now =0
t=1; %iteration counter
fprintf(' cacluating Non-Negative EVecs for alpha = %d and beta = %d\n',alpha1,beta);
while (1)

    for r=1:k
        
        for s=1:d
            c1 = C(s,:)*u(:,r) -  C(s,s)*u(s,r) - alpha1*( ( u(s,1:k) - u(s,r) )*( u(:,r)'*u(:,1:k) - u(s,r)*u(s,1:k) )' ) - beta;
            
            c2 = C(s,s)+ alpha1 - alpha1*(u(:,r)'*u(:,r)-u(s,r)^2) - alpha1*(u(s,1:k)*u(s,1:k)' -u(s,r)^2) ;
            
            poly = [-alpha1 0 c2 c1];
            roots_df = roots(poly);
            roots_df ( imag(roots_df)~=0 ) = 0 ; % only real!
            roots_df ( roots_df<0 ) = 0 ; % only positive
            roots_check=[0;roots_df];
            objective = -alpha1/4*roots_check.^4 + c2/2*roots_check.^2 + c1*roots_check;
            [~,I] = max(objective);

            u(s,r) =  roots_check(I);
            if ( u(s,r)>1e12)
                    fprintf('error! to big!')
                    break;
            end
        end  %end of for. claculated all the column  vector
        
    end   %end of 2nd for. clculated the entire row
    
         t=t+1;  %t is the time counter. for exiting the while (iterations)
         if mod(t,100)==0
             fprintf('iteration #%d \n',t)
         end
         if (t>2e3) || ( norm(u_prev-u,'fro')<0.0001 )
                 fprintf('Iterations at %d Breaking out...\n',t-1)
             break; %get out of the while if iterations too big or convergance satisfied
         end
        u_prev=u;
% %         ------u =  aux;
  
end %end of while
  fprintf('Finished calcs, printing new EigenVectors\n');
%    figure; plot(u(:,1:8))
%    for q=1:k
%       figure;
%       plot( u(:,q))
%       title(['Evec with alpha=',num2str(alpha1)])
%    end
% hold on
% plot( u(:,2),'r')
%  hold on
% plot( u(:,5),'g')   
end