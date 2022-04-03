function [p] = half_hermite_shift_poly(x, n, u)
% hermite_poly -- Evaluates Hermite polynomials
%
% p = hermite_poly(x, n, [[hmean=0, hvar=1, normalize=0]])
%
%     Evaluates the degree-n Hermite polynomials as the location x. This code is
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

[alpha, beta] = half_hermite_shift_recurrence_symb(N+2,u);

% Begin recurrence
p = zeros([length(x) length(n)]);
p0 = zeros([length(x) 1]);
p1 = zeros([length(x) 1]);

temp = (1 + erf(u))*sqrt(pi)/2;
% temp_s = sym(temp);
% temp_s = temp_s * sqrt(pi_s)/2;

p0(:) = ones(xsize)/sqrt(temp);
p1(:) = (x-alpha(1)).*p0(:);
p1(:) = p1(:)/sqrt(beta(2));
% p1(:) = (x-1/sqrt(pi))/sqrt(sqrt(pi)/4-1/2/sqrt(pi));

alpha = alpha(2:end);
beta = beta(2:end);

p(:,n==0) = p0;
p(:,n==1) = p1;

% Loops over the recurrence
for q = 2:N
  temp = p1;
  p1 = 1/sqrt(beta(q))*((x-alpha(q-1)).*p1 - sqrt(beta(q-1))*p0);
  p0 = temp;

  p(:,n==q) = p1;
end

% Format output
if length(x)==1
  p = reshape(p, nsize);
elseif length(n)==1
  p = reshape(p, xsize);
end
return