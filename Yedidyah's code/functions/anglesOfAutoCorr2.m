
function  [angles,minAngles,secAngle,minDist] =anglesOfAutoCorr2(xCoor,yCoor,xx,cross_map, paint)
centers= [(numel(xx)+1)/numel(xx)*max(xx), (numel(xx)+1)/numel(xx)*max(xx)];
clear coords
coords (1,:) = xCoor';
coords (2,:) = yCoor';

%take only right half hemsphire
xCoords = coords(1,:);
limit_x = 0;
flag=1;
while (flag==1)
    rightSide_inds = find((xCoords>=centers(1))- limit_x );
    rightSide_coords = coords(:,rightSide_inds);
    
    Dists = dist(centers,coords(:,rightSide_inds)) ; %#ok<FNDSB>
    
    [~, indDists] = sort(Dists);
    if (Dists(indDists(1)) < 0.5) %if the center is the first one
        first_ind = 2;
    else
        first_ind = 1;
    end
    
    minDist = Dists(indDists(first_ind)); % take the 2nd from the start, since the 1st one is always zero.
    
    maxM = min(first_ind+5,length(indDists)); %choose the 1st 6 (if exists) distances
    top6Dists = rightSide_coords(:,indDists(first_ind:maxM));
    if  length(top6Dists)<4
        limit_x = limit_x +5;
    else
        flag = 0;
    end
end

 MaxN = min(length(top6Dists),first_ind+2);
% [~, indDists3top] = sort(top6Dists(2,:));
% top3Dists = top6Dists(:,indDists3top(1:MaxN));
top3Dists = rightSide_coords(:,indDists(first_ind:MaxN));

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
%Check if angles "make sense". sometimes this can lead to error
% if minAngles>35 
%     minAngles = 1000; %ignore these...!
% end

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
