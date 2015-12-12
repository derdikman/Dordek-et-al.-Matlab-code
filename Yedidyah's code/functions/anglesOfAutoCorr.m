
function  [angles,minAngles,secAngle,minDist] =anglesOfAutoCorr(xCoor,yCoor,xx,cross_map, paint)
centers= [(numel(xx)+1)/numel(xx)*max(xx), (numel(xx)+1)/numel(xx)*max(xx)];
clear coords
coords (1,:) = xCoor';
coords (2,:) = yCoor';

Dists = dist(centers,coords) ;

[~, indDists] = sort(Dists);
if (Dists(indDists(1)) < 0.5) %if the center is the first one
    first_ind = 2;
else
     first_ind = 1;
end
    
minDist = Dists(indDists(first_ind)); % take the 2nd from the start, since the 1st one is always zero.

maxM = min(max(indDists),5);
top6Dists = coords(:,indDists(first_ind:maxM));

MaxN = min(length(top6Dists),3);
[~, indDists3top] = sort(top6Dists(2,:));
top3Dists = top6Dists(:,indDists3top(1:MaxN));

if paint ==1
    figure; imagesc(2*xx,2*xx , cross_map), colormap jet, hold on
    for i=1:length(top3Dists)
        plot(top3Dists(1,i),top3Dists(2,i),'X');
        hold on
        plot([centers(1),top3Dists(1,i)],[centers(2),top3Dists(2,i)])
        hold on
    end
end

a1 = centers(1) -  top3Dists(1,:);
a2 = centers(2) -  top3Dists(2,:);
angles = abs( atand(a2./a1))';
 [ minVals, idd] = min(min(abs(angles - 90), angles));
% minAngles = angles(idd);

minAngles = min(min(angles), min(abs(angles-90)));
% if (sum (angles<0) > 1 )
%     angles = atand(a2./a1) +180;
% end

%2nd angle
idx = 1:length(angles);
sortedAngles= angles(idx~=idd);
secAngle  = max(90-sortedAngles(1) ,sortedAngles(1));


% angles_raw = atand(a2./a1);
% sortedAngles = sort(mod(angles_raw-89.999999,90));
% secAngle = max(sortedAngles);
end
