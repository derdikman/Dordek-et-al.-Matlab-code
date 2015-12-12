function  [AHex, ASq]=ODEHornikCalc(Mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Hornik's ODE, when given a covarince matrix. Does the outputs come out Hexagonal when initilized randomly?

Mat = r_s;

NmEC = 100;
L=length(Mat);
meanInputMat = mean(Mat(:,1:L),2);
InputMatFixed = Mat(:,1:L) - repmat(meanInputMat,1,L);
C = InputMatFixed*InputMatFixed';
Cov1 = (1/L)*C;
N = sqrt(length(Cov1));
figure; imagesc(Cov1), colormap jet
for i=1:2
    % initilize A
    i = 2 %%%%%%%%%%%%%%%%%%%%%%%CHECK 
    
    if (i==1) %%%%%%%!!!!!!!FIX THE SIZE! %%% TO DO!!! (needs to be 100X400 here
        A =-.1+ .2*rand(NmEC,length(Cov1) ); %random values - negative as well
    else
        A = rand( NmEC,length(Cov1)  ) ; %only positive
    end
    normFactor = repmat(sum(A,2),1,size(Cov1,2)); % normelize A
    A = A./normFactor;
    W  = zeros(NmEC);
    %time index
    tind =1;
    simdur = 5e4;
    %figure
    %loop
    while (tind<simdur)
        %the ODE
%         Q = inv(eye(size(W) ) -W) ;
%         R = Q*A*Cov1*A'*Q';
%         ADot = (Q)*A*Cov1 - eye(NmEC).*(R)*A; %
        
        R = A*Cov1*A';
         ADot = A*Cov1 - eye(NmEC).*(R)*A;
        % set epsilon to fit the conditions HORNIK states, see condition A2, page 233 in his paper from 1991
        epsilon = 1./(tind*.3+1e4);
        %update A
        A = A + epsilon*ADot;
        %cut A if we want HEXAGONAL results
        if i==2
            A(A<0)=0;
        else
            A(A<-10) = -10;
            A(A>10) = 10;
        end
        %update W
%         WDot =tril( (ones(size(W)) - eye(size(W))).*R,-1);
%         W = W - epsilon*WDot;
%          W(W<-.1) = -.1;
%          W(W>.1) = .1;
        
        tind  = tind +1;
        tind
        
%         updateRate(tind) = sum(sum( epsilon*ADot))/(size(Cov1,2))^2;
            if (mod(tind,100)==0)
        
                pause (0.0001)
                imagesc(reshape(A(1,:),N,N)), colormap jet
                drawnow
            end
    end
    if i==2
        AHex = A;
    else
        ASq = A;
    end
end
end
% 
% figure; imagesc(reshape(A(5,:),N1,N1)),colormap jet
%A1 = reshape(A(10,:),20,20);
% n_mat=mulMatCross(zeros(20),zeros (20));
% cross_mapA1=Cross_Correlation(A1,A1,n_mat);
% figure;
% imagesc(cross_mapA1), colormap jet
% Gridness60 = findGrindess( cross_mapA1);
% Gridness90 = findGrindess90( cross_mapA1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear pc
counter = 1;
%create aux vec
placeCentersAux = placeCenters1;
%start with the frist one
 pc(:,counter) =placeCentersAux(:,1);
for i =1: length(placeCenters1)
      %init distance min
    dmin = 100;
    for j=1: length(placeCenters1)
        %this indicates that this place cell is already taekn
        if placeCentersAux(1,j) == 1000
            continue;
        end
        %with modulo diff
        diffX = min(mod(pc(1,counter)  - placeCentersAux(1,j),maxx),mod( placeCentersAux(1,j)-pc(1,counter) ,maxx) );  
       diffY = min(mod(pc(2,counter)  - placeCentersAux(2,j),maxx),mod( placeCentersAux(2,j)-pc(2,counter) ,maxx) );  
            
        dTmp = diffX.^2  +diffY.^2 ;
       if dTmp<dmin && dTmp>0
           dmin  = dTmp ;
           indMinDist = j;
       end
    end
    %dMinDist(i,indMinDist) = dmin;
    %save the indeces
    IDX(i) = indMinDist;
    pc(:,counter+1) = placeCentersAux(:, indMinDist);
    %make sure we wont touch this place again
    placeCentersAux(:, indMinDist) = [ 1000; 1000];
    counter = counter+1;
end

