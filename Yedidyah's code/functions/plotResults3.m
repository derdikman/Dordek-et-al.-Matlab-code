function [  ] = plotResults3( unitMaps,n_mat,m,title1)

if (size(n_mat,1)~=size(unitMaps,1)*2-1)
    fprintf('n_mat not the right size. recalculating..\n')
    sz = size(unitMaps,1) ;
    n_mat = mulMatCross(zeros(sz,sz), zeros (sz,sz) ) ;
end

e=1;
for p=1:m
    figure;
    for q=1:9
        subplot(3,3,q)
        imagesc(unitMaps(:,:,e));
        title(['#PC ' ,num2str(e)])
        axis equal
        e=e+1;
    end
    suptitle(['Averaged activity of ',title1])
end
%return;
e=1;
str = 'Autocorrelation of PC#';
for p=1:m
    figure;
    for q=1:9
        subplot(3,3,q)
        mat=unitMaps(:,:,e) ;
        mat_cross=Cross_Correlation(mat,mat,n_mat);
        imagesc(mat_cross)
        title( sprintf( '%s %d', str, e ) );
        axis equal
        e=e+1;
    end
    suptitle(['Autocorrelation of ',title1])
end


end

