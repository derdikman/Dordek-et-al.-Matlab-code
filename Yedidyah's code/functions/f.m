function [ x_out] = f( x,n)
x(x<0)=0;
x_out = g( x,n);
end