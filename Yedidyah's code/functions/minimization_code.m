


%%
%init w0
CStruct = load('Cov.mat');
Cov = CStruct.Cov;
 M = 196;
 Cov = Cov(1:M,1:M) +abs(min(min(Cov(1:M,1:M)))); %dosnt work on negative??
% %Cov = rand(100);
% % p = randperm(size(C,1)); %create randpermutation of the columns of Cov
% % C =Cov(p,:);

C = Cov;
figure; imagesc(C), colorbar
alpha1 = 1e1;
d = size(C,1) ; %size of each vector
k=3; %number of rows (vectors)
u = randn(d,k);          %initial guess for u!!!
beta = 0; %for now =0
t=1; %iteration counter
%clear aux
% aux(:,1) = rand(1,size(Cov,1)); %initial guess for aux
% aux(2,:) = rand(1,size(Cov,1)); %initial guess for aux

while (1)

    for r=1:k
        
        for s=1:d
            c1 = C(s,:)*u(:,r) -  C(s,s)*u(s,r) - alpha1*( ( u(s,1:k) - u(s,r) )*( u(:,r)'*u(:,1:k) - u(s,r)*u(s,1:k) )' ) - beta;
            
            c2 = C(s,s)+ alpha1 - alpha1*(u(:,r)'*u(:,r)-u(s,r)^2) - alpha1*(u(s,1:k)*u(s,1:k)' -u(s,r)^2) ;
            
            poly = [-alpha1 0 c2 c1];
            roots_df = roots(poly);
            roots_df ( imag(roots_df)~=0 ) = 0 ; % only real!
            nonnegroot = max(roots_df,0 );% only positive!
            aux(s,r) =    max(max(-alpha1/4*nonnegroot.^4+c2/2*nonnegroot.^2+c1*nonnegroot),0);
            
    
            
        end  %end of for. claculated all the column  vector
        
    end   %end of 2nd for. clculated the entire row
    
         t=t+1  %t is the time counter. for exiting the while (iterations)
         if (t>1e2) || ( abs( norm(aux,'fro') - norm(u,'fro'))<0.001)
            break; %get out of the while if iterations too big or convergance satisfied
         end
        
        u =  aux;
    
end %end of while

figure; 
plot( u(:,2),'-.')
        
        %% for K=1 /only 1 EigVec
    
        CStruct = load('Cov.mat');
Cov = CStruct.Cov;

M = 196;
C = Cov(1:M,1:M) +abs(min(min(Cov(1:M,1:M)))); %dosnt work on negative??
%Cov = rand(100);
% p = randperm(size(C,1)); %create randpermutation of the columns of Cov
% C =Cov(p,:);

%C = Cov;
            w = rand(2,size(C,1));
        t=2;
        d=M;
        alpha1 = 1e2;
        
        while  (t<1e3) && (norm( w(t,:) - w(t-1,:))>0.001)
            for s=1:d
                c1 = C(s,:)*w(t,:)' - C(s,s)*w(t,s)  ;
                c2 = C(s,s)+ alpha1 - alpha1*(w(t,:)*w(t,:)'-w(t,s)^2);
                poly = [-alpha1 0 c2 c1];
                roots_df = roots(poly);
                nonnegroot = roots_df(roots_df>0);
                w(t+1,s)=    max(max(-alpha1/4*nonnegroot.^4+c2/2*nonnegroot.^2+c1*nonnegroot),0);
                %         nonnegroot= abs(max(vertcat(roots_df,0))); %what happens it the root is imaginary??
                %         nonnegroot = roots_df;
                %         w(t+1,s)= -alpha1/4*nonnegroot.^4+c2/2*nonnegroot.^2+c1*nonnegroot;
                if ( w(t+1,s)>1e12)
                    fprintf('error! to big!')
                    break;
                end
            end
            
            t=t+1
        end
        figure; plot(  w(t,:),'-.')
        figure; plot(w(t,:),'.')
hold on
plot(w(t-3,:),'r.')
hold on
 plot(w(t-2,:),'g')
