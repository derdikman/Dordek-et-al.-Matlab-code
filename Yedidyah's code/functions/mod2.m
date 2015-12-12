function val = mod2(X,Y)
% this function returns regular modulo BUT if the value returned is 0 it
% goes back to the boundry value (for IMAGESC indeces cacls)
val = mod(X,Y);
if val == 0
    val = Y;
end

end