function [outputMap, autocorrMap, Gridness60, Gridness90,minAngles,minDists] = Maps(Struct,weights)

%unpack vars
placeCenters = Struct.placeCenters;
NN=Struct.Resolution; %Map Resolution
psiSat =Struct.psiSat;

%limit size of Maps...
if size(weights,2)>100
    k = 100;
else 
    k = size(weights,2);
end

G= zeros(NN,NN,k);
n_mat=mulMatCross(zeros(size(G,1)),zeros (size(G,1)));
xx=linspace(0,Struct.maxx,NN);
yy=linspace(0,Struct.maxy,NN);
[R,C1] = ndgrid(xx, yy);

sigma_x = Struct.PcSize;
sigma_y = Struct.PcSize;
theta = 2*pi*rand(1,length(placeCenters)); %Gives a "tilt" to the place cell... not used when sigma_x =sigma_y
a_std = cos(theta).^2/2/sigma_x.^2 + sin(theta).^2/2/sigma_y.^2;
b_std = -sin(2*theta)/4/sigma_x.^2 + sin(2*theta)/4/sigma_y.^2 ;
c_std = sin(theta).^2/2/sigma_x.^2 + cos(theta).^2/2/sigma_y.^2;
% place cells sombrero hat - 2nd Gaussians
sigma_x2 = 2*sigma_x;
sigma_y2 =2*sigma_y;
theta2 = 2*pi*rand(1,length(placeCenters));
a_std2 = cos(theta2).^2/2/sigma_x2.^2 + sin(theta2).^2/2/sigma_y2.^2;
b_std2 = -sin(2*theta2)/4/sigma_x2.^2 + sin(2*theta2)/4/sigma_y2.^2 ;
c_std2 = sin(theta2).^2/2/sigma_x2.^2 + cos(theta2).^2/2/sigma_y2.^2;

fprintf(['Calculating maps, input: ',Struct.placeCellType,'\n'])
clear cross_map
for q=1:k
    for i=1:size(placeCenters,2)
        diff_x=min(mod(R-placeCenters(1,i),Struct.maxx),mod(placeCenters(1,i)-R,Struct.maxx));
        diff_y=min(mod(C1-placeCenters(2,i),Struct.maxy),mod(placeCenters(2,i)-C1,Struct.maxy));
        % Place cell activation at this location
        squareDists = (a_std(i).*(diff_x).^2 +2*b_std(i).*(diff_x).*(diff_y) + c_std(i).*(diff_y).^2 );
        squareDists2 = (a_std2(i).*(diff_x).^2 +2*b_std2(i).*(diff_x).*(diff_y) + c_std2(i).*(diff_y).^2 );
        
        A = sigma_x^2/sigma_x2^2;
        if strcmp(Struct.placeCellType,'DOG')
            val(:,:,i) =  psiSat*(exp(-squareDists) - A*exp(- squareDists2));
        elseif  strcmp(Struct.placeCellType,'Gaussian')
            val(:,:,i) =  psiSat*(exp(-squareDists));
        elseif  strcmp(Struct.placeCellType,'Disk')
            val(:,:,i) = Struct.DisksData(:,:,i);
        end
        G(:,:,q)= G(:,:,q)+  weights(i,q) *val(:,:,i);
    end
    
    cross_map(:,:,q) = Cross_Correlation( G(:,:,q), G(:,:,q),n_mat);
    Gridness60(q) = findGrindess( cross_map(:,:,q));
    Gridness90(q) = findGrindess90( cross_map(:,:,q));
    
    fprintf(['Calculating output #',num2str(q),'...\n'])
end

outputMap = G;
autocorrMap = cross_map;

Map = cross_map;
paint = 0;
thresh = 0;
rad = 0.75*NN/Struct.maxx;
angles33 = zeros(3, size(Map,3));
field1 = 'centers';
field2 = 'distances';
H = fspecial('disk',10);
fprintf('Calculating Grid-spacing and alignment \n')
for qq = 1:size(Map,3)
    [localMaxIndF] = transpose ( localMaxFinder(Map(:,:,qq),thresh,rad,H,paint));
    cent(qq).field1 = localMaxIndF(:);
    xCoor =  cent(qq).field1(1:2:end).*Struct.maxx/(NN);
    yCoor =  cent(qq).field1(2:2:end).*Struct.maxy/(NN);
    [angles33(:,qq),minAngles(qq),secAngle(qq),minDists(qq)] =anglesOfAutoCorr2(xCoor,yCoor,xx,Map(:,:,qq), paint);
end





end