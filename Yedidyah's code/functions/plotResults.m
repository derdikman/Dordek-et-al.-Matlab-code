function [  ] = plotResults( totalActivity,n_mat,m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%totalActivity(totalActivity<0)=0;
str = 'Autocorr'' #';

if size(totalActivity,3) <9
    for i=1:size(totalActivity,3)
        figure;
               imagesc(totalActivity(:,:,i)./(eps)); colormap jet
        title([' output  #' ,num2str(i)])
              figure;      
              mat=totalActivity(:,:,i)./(eps) ;
        mat_cross=Cross_Correlation(mat,mat,n_mat);
        imagesc(mat_cross), colormap jet
       title( sprintf( '%s %d', str, i ) );
        
        
        
    end
else

e=1;
for p=1:m
    figure;
   ha = tight_subplot(2,4,[.03 .03],[.1 .1],[.01 .01]);
    for q=1:8
        axes(ha(q))
      %  subplot(3,3,q)
        imagesc(totalActivity(:,:,e)./(eps)); colormap jet
        title([' output  #' ,num2str(e)])
     
    % axis auto
        axis off
        axis equal
        e=e+1;
    end
   % suptitle(['Averaged activity of ',title1])
    % title2 = ({'Averaged activity of neural networks'' outputs';'Weights are constrained'});
%suptitle( (title2))
end

e=1;%Autocorrelation of 
for p=1:m
    figure;
       ha = tight_subplot(2,4,[.03 .03],[.1 .1],[.01 .01]);

    for q=1:8
                axes(ha(q))
%        subplot(3,3,q)
        mat=totalActivity(:,:,e)./(eps) ;
        mat_cross=Cross_Correlation(mat,mat,n_mat);
        imagesc(mat_cross), colormap jet
       title( sprintf( '%s %d', str, e ) );
       axis equal
        axis off
        e=e+1;
    end
      %   title2 = ({'Autocorrelation of neural networks'' outputs';'Weights are constrained'});
%suptitle( (title2))
%      title2 = ({'Averaged activity of neural networks'' outputs';''});
% suptitle( strcat(title2,title1))
   % suptitle(['Autocorrelation of ',title1])
end

end
end
% gap = round(simdur/300);
% figure
%            ha = tight_subplot(3,3,[.04 .04],[.1 .1],[.01 .01])
% 
% for a=1:9
%              axes(ha(a))
%     %subplot(3,3,a)
%     numc = round((a-1)*4 +2);
%     imagesc(G1(:,:,numc))
%    title(['After ',num2str(round(numc*gap/1e3)),'e3 time steps'])
%     axis equal
%        axis off
% end
%  suptitle({'Evolution in time of neural network output';'Weights are constrained'})
% 
% 
% gap = round(simdur/300);
% figure
%            ha = tight_subplot(3,3,[.04 .04],[.1 .1],[.01 .01])
% 
% for a=1:9
%              axes(ha(a))
%     %subplot(3,3,a)
%     numc = round((a-1)*5 +2);
%             mat_cross=Cross_Correlation(G1(:,:,numc),G1(:,:,numc),n_mat);
%         imagesc(mat_cross)
%    title(['After ',num2str(round(numc*gap/1e3)),'e3 time steps'])
%     axis equal
%        axis off
% end
%  suptitle({'Evolution in time of the output''s autocorrelation ';'Weights are constrained'})



%% BACKUP
% 
% 
% for p=1:m
% figure;
% 
% subplot(331)
% imagesc(totalActivity(:,:,q+1)./(eps+occupancyTime));
% %title([title2 ,num2str(q+1)])
% axis equal
% subplot(332)
% imagesc(totalActivity(:,:,q+2) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+2)])
% axis equal
% subplot(333)
% imagesc(totalActivity(:,:,q+3) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+3)])
% axis off
% subplot(334)
% imagesc(totalActivity(:,:,q+4) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+4)])
% axis off
% subplot(335)
% imagesc(totalActivity(:,:,q+5) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+5)])
% axis off
% subplot(336)
% imagesc(totalActivity(:,:,q+6) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+6)])
% axis off
% subplot(337)
% imagesc(totalActivity(:,:,q+7) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+7)])
% axis off
% subplot(338)
% imagesc(totalActivity(:,:,q+8) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+8)])
% axis off
% subplot(339)
% imagesc(totalActivity(:,:,q+9) ./ (eps + occupancyTime))
% %title([title2 ,num2str(q+9)])
% %axis equal
% set(gca,'position',[0 0 1 1],'units','normalized')
% axis off
% title2 = ({'Averaged activity of ';''});
% suptitle( strcat(title2,title1))
% q=q+9;
% end
% e=1;
% q=0;
% %str = 'Autocorrelation ';
% str = '';
% 
% 
% for p=1:m
% figure;
% subplot(331)
% mat(:,:,e)=totalActivity(:,:,q+1) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+1;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% 
% e=e+1;
% subplot(332)
% mat(:,:,e)=totalActivity(:,:,q+2) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+2;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% 
% e=e+1;
% subplot(333)
% mat(:,:,e)=totalActivity(:,:,q+3) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+3;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(334)
% mat(:,:,e)=totalActivity(:,:,q+4) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+4;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(335)
% mat(:,:,e)=totalActivity(:,:,q+5) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+5;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(336)
% mat(:,:,e)=totalActivity(:,:,q+6) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+6;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(337)
% mat(:,:,e)=totalActivity(:,:,q+7) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+7;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(338)
% mat(:,:,e)=totalActivity(:,:,q+8) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+8;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% e=e+1;
% 
% subplot(339)
% mat(:,:,e)=totalActivity(:,:,q+9) ./ (eps + occupancyTime);
% mat_cross(:,:,e)=Cross_Correlation(mat(:,:,e),mat(:,:,e),n_mat);
% imagesc(mat_cross(:,:,e))
% num=q+9;
% %title( sprintf( '%s%s %d', str, title2, num ) );
% axis equal
% title2 = ({'Autocorrelation of  ';''});
% suptitle( strcat(title2,title1))
% e=e+1;
% 
% q=q+9;
% end