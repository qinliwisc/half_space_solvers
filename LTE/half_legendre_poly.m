function [p] = half_legendre_poly(x, n)
% half_lagendre_poly -- Evaluates Legendre polynomials
%
%
%     Evaluates the degree-n Legendre polynomials as the grid points x. This code is
%     `vectorized' in both x and n:
%
%        hermite_poly(x, n) if x is an array and n is a scalar returns the
%            degree-n Hermite polynomial evaluated at x, and the output array is
%            of the same size as x
%

% Pre-processing, resizing:
xsize = size(x); x = x(:); 
nsize = size(n); n = n(:);
N = max(n);

[alpha, beta] = half_legendre_recurrence(N);

% Begin recurrence
p = zeros([length(x) length(n)]);
p0 = zeros([length(x) 1]);
p1 = zeros([length(x) 1]);

p0(:) = ones(xsize);
p1(:) = (2*x-1)*sqrt(3);

p(:,n==0) = p0;
p(:,n==1) = p1;

% Loops over the recurrence
for q = 2:N
  temp = p1;
  p1 = ((x-alpha(q-1)).*p1 - beta(q-1)*p0)/beta(q);
  p0 = temp;

  p(:,n == q) = p1;
end

% Format output
if length(x)==1
  p = reshape(p, nsize);
elseif length(n)==1
  p = reshape(p, xsize);
end
return