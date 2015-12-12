function [c,ceq] = mycon(x)
c =[];   % Compute nonlinear inequalities at x.
ceq = norm(x)-1;   % Compute nonlinear equalities at x.
end