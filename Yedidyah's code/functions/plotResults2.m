function [  ] = plotResults2( unitMaps,n_mat,m,title1)

if (size(n_mat,1)~=size(unitMaps,1)*2-1)
    fprintf('n_mat not the right size. recalculating..\n')
    sz = size(unitMaps,1) ;
    n_mat = mulMatCross(zeros(sz,sz), zeros (sz,sz) ) ;
end

q=0;
for p=1:m
figure;
subplot(331)
imagesc(unitMaps(:,:,q+1));
title([title1 ,num2str(q+1)])
axis equal
subplot(332)
imagesc(unitMaps(:,:,q+2) )
title([title1 ,num2str(q+2)])
axis equal
subplot(333)
imagesc(unitMaps(:,:,q+3) )
title([title1 ,num2str(q+3)])
axis equal
subplot(334)
imagesc(unitMaps(:,:,q+4) )
title([title1 ,num2str(q+4)])
axis equal
subplot(335)
imagesc(unitMaps(:,:,q+5) )
title([title1 ,num2str(q+5)])
axis equal
subplot(336)
imagesc(unitMaps(:,:,q+6) )
title([title1 ,num2str(q+6)])
axis equal
subplot(337)
imagesc(unitMaps(:,:,q+7) )
title([title1 ,num2str(q+7)])
axis equal
subplot(338)
imagesc(unitMaps(:,:,q+8) )
title([title1 ,num2str(q+8)])
axis equal
subplot(339)
imagesc(unitMaps(:,:,q+9) )
title([title1 ,num2str(q+9)])
axis equal
suptitle(['Averaged activity of ',title1])
q=q+9;
end
e=1;
q=0;
str = 'Autocorrelation ';
for p=1:m
figure;
subplot(331)
mat(:,:,e)=unitMaps(:,:,q+1) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+1;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal

e=e+1;
subplot(332)
mat(:,:,e)=unitMaps(:,:,q+2) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+2;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal

e=e+1;
subplot(333)
mat(:,:,e)=unitMaps(:,:,q+3) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+3;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(334)
mat(:,:,e)=unitMaps(:,:,q+4) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+4;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(335)
mat(:,:,e)=unitMaps(:,:,q+5) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+5;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(336)
mat(:,:,e)=unitMaps(:,:,q+6) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+6;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(337)
mat(:,:,e)=unitMaps(:,:,q+7) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+7;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(338)
mat(:,:,e)=unitMaps(:,:,q+8) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+8;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
e=e+1;

subplot(339)
mat(:,:,e)=unitMaps(:,:,q+9) ;
mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
imagesc(mat_cross(:,:,e))
num=q+9;
title( sprintf( '%s%s %d', str, title1, num ) );
axis equal
suptitle(['Autocorrelation of ',title1])
e=e+1;

q=q+9;
end

end

