function gamma = bestGamma(r1,r2,xyc,gamma1,NN,X1,Y1)
I=zeros(NN);
gamma2F = r1^2/(r2^2-r1^2);
maxx = max(max(X1));
maxy = max(max(Y1));
diff_x=min(mod(X1-xyc(1),maxx),mod(xyc(1)-X1,maxx));
diff_y=min(mod(Y1-xyc(2),maxx),mod(xyc(2)-Y1,maxy));
circleEq =diff_x.^2+ diff_y.^2;
I  (circleEq<r1^2  )  = gamma1;
I  ( (circleEq>=r1^2  )  &  (circleEq<=r2^2  ))  = - gamma2F;
gamma = sum(sum(I==gamma1))/sum(sum(I==-gamma2F));

end

%% Second version (not modulo
% I2=zeros(NN);
% gamma2F = r1^2/(r2^2-r1^2);
% maxx = max(max(X1));
% maxy = max(max(Y1));
% I2  ( ((X1-xyc(1)).^2 + (Y1-xyc(2)).^2) <  r1^2 )  = gamma1;
% I2( ( ((X1-xyc(1)).^2 + (Y1-xyc(2)).^2) >=  r1^2 ) &...
%     (((X1-xyc(1)).^2 + (Y1-xyc(2)).^2) <= r2^2) ) =-gamma;
% gamma = sum(sum(I2==gamma1))/sum(sum(I2==-gamma2F));

