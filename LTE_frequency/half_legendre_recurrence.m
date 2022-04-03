function [alpha,beta] = half_legendre_recurrence(N)
% recurrence -- recurrence coefficients for Legendre polynomials
%
% [a,b] = recurrence(N)
%     Calculates the first N recurrence coefficients for the hermite
%     polynomials.  


if nargin == 0
    N = 15;
end

alpha = 0.5*ones(N+1,1);
beta = [0:1:N];
beta = (beta+1)./sqrt(2*beta+1)./sqrt(2*beta+3)/2;

return