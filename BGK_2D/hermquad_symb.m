function [X,W] = hermquad_symb(n) % the method presented in Golub's paper
if nargin==0
    n = 15;
end

[alpha] = hermite_recurrence_symb(n);

A = diag(alpha(2:end), 1) + diag(alpha(2:end), -1);

[V,X]=eig(A);
X = diag(X);
W = vpa(V(1, :)'.^2)*sqrt(pi);

X = double(X);
W = double(W);