function [X,W]=half_legendre_quad(n) % the method presented in Golub's paper
if nargin==0
    n = 15;
end

[alpha,beta] = half_legendre_recurrence(n);

A = diag(alpha(1:end)) + diag(beta(1:end-1), 1) + ...
    diag(beta(1:end-1), -1);

[V,X]=eig(A);
X = diag(X);
W = vpa(V(1, :)'.^2);

W = double(W);

return