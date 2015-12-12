

%[x,fval,exitflag,output] =  fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

%C= randn(100);
%Cov = load('Cov.mat');
% C = Cov.Cov(:,randperm(625));
M=180;
C = Cov(1:M,1:M) +abs(min(min(Cov(1:M,1:M))));

p = randperm(size(C,1)); %create randpermutation of the columns of Cov
C =C(p,:);

clear x5
clear x0
for e=1:1
x0(:,e)=rand(size(C,1),1);
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunEvals',100000);
lb =zeros(size(C,1),1); % -Inf(size(C,1),1);
ub = 1e5*ones(size(C,1),1);%Inf(size(C,1),1);
[x5(:,e) ,fval,exitflag,output] = fmincon(@(x)myfun(x,C),x0(:,e),[],[],[],[],lb,ub,@(x)mycon(x),options);
e
end
 figure; plot(x5)
 figure; plot(x0)
 
figure; plot(x5(:,4),'b.'), hold on,plot(x0(:,1),'r.')


%% QUADPROG
options = optimoptions(@quadprog,'MaxIter',10000);
[x6,fval,exitflag,output] =quadprog(-1*C,[],[],[],[],[],lb,ub,x0,options);
