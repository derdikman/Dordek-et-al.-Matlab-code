function [next_x, next_y,curDir] = movement4PCA(x,y,curDir,maxx,dirStd,velocity)

        directions = [0,0;
                               0,1;
                               0,-1;
                               1,0
                               0,-1;
                               1,1;
                               -1,-1;
                               -1,1;
                               1,-1];
         randDir = datasample(1:9,1);   
      %  randDir = round(1+ mod(curDir + dirStd*rand,8));
        next_x = mod(x + velocity*directions(randDir,1),maxx);
         next_y =mod( y + velocity*directions(randDir,2),maxx);
        curDir = randDir;
        
end