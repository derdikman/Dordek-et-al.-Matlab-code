function [localMaxIndF] = localMaxFinder(Image,thresh,rad,H, paint)


AA = Image;
localMaxLen =0;
localMaxInd = [1000 1000];
[M, N] = size(AA);
 %clear localMaxInd

 

 %AA = imfilter(Image,H,'replicate');
  if (paint ==1)
 figure;
  imagesc(AA);colormap jet
% %axis image
 hold on;
 end
%search for local maximas. Don't start from corner to avoid effects if not
%in modulo. later?
for i = 1:size(AA,1)-1
    for j = 1:size(AA,2)-1
              
       if AA(i,j) > AA(mod2(i+1,M),j) && ...
           AA(i,j) > AA(mod2(i-1,M),j) && ...
                 AA(i,j) > AA(i,mod2(j+1,N)) && ...
                 AA(i,j) > AA(i,mod2(j-1,N)) && ...
                   AA(i,j) > AA(i+1,mod2(j+1,N)) && ...
                     AA(i,j) > AA(i+1,mod2(j-1,N)) && ...
                          AA(i,j) > AA(mod2(i-1,N),j+1) && ...
                     AA(i,j) > AA(mod2(i-1,N),mod2(j-1,N)) && ...
                  AA(i,j) >=thresh && ... %neutrelize negative vals/under threshold
                  sum( ( localMaxInd(:,1) - i).^2 + (localMaxInd(:,2) - j).^2  > rad^2) >= (localMaxLen)  %checkiing that the new position isnt too close. avoid error
                %(sum(localMaxLen) >1) && 
              if (paint ==1)
              plot(j,i,'x');
              end
              localMaxInd(localMaxLen+1,1) = i;
              localMaxInd(localMaxLen+1,2) = j;
              localMaxLen = localMaxLen+1;    

       end
    end
end
% if localMaxLen==0  %if there's non
%     break; 
% end
localMaxInd(:,[1 2]) = localMaxInd(:,[2 1]); 
localMaxVal = diag(AA(localMaxInd(:,2),localMaxInd(:,1))); %get the values of the local maximas
[val, ind]= sort(localMaxVal,'descend');
    localMaxIndF1 = localMaxInd(ind,:);
    localMaxValF1=val;
   %%   %look for unique and distanct locations only!
     for i=1:length(val)
         for j=1:length(val)
          p(i,j) = pdist( [localMaxIndF1(i,:); localMaxIndF1(j,:)]);
         end
     end
     
     flag =100* ones(length(val));
     
     for i=1:length(val) %row
         for j=1:length(val) %column
             if (j==i)  %don't check self distnace
                 continue;
             end
             if (p(i,j) < 5)   %radius of "too" close maximas (not too local..) - Can change when resolution changes
                  if sum(sum( ismember(flag,p(i,j)) )) == 0  %check if the flag was raised for this maxima or its neighbor
                        flag(i,j) = p(i,j);  %if not - place it in flag for next loop to skip it.
                 else
                     break; %if so - skip to the next maxima and don't take this one
                 end
             end
                         
             end
             if (j==length(val))   %if we managed to go over all distances from maxima without being kicked out - OK!
                indFinal (i) = i;
             end
         end
      localMaxIndF = localMaxIndF1(indFinal~=0,:); %take only the unique/with no neigbors ones...
      localMaxValF = val(indFinal~=0);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% get from the distances from each point to its neighbor - devide by 2
%     % and round for finding uniqeness
%     pF = round(p(:,1)./2);
%     [~ ,indUniq]= unique( pF);
%     %%
%      localMaxIndF =localMaxIndF1(indUniq,:);
%      localMaxValF = val(indUniq);
end