function [f] = myfun(x,C)
f = -x'*C*x;            % Compute function value at x
%g=-2*C*x;
end