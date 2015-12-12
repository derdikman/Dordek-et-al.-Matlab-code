function [ x_out] = g( x,n)
if norm(x,2)==0
    x_out=zeros(length(x),1);
else
    x_out = sqrt(n)*x/(norm(x,2)+eps);
end
end