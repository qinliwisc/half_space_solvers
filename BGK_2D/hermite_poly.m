function [p] = hermite_poly(x, N)

if nargin == 0
    x = [-1:0.01:1]; N = 5;
end

% Pre-processing, resizing:
xsize = size(x); x = x(:);

% Begin recurrence
p = zeros([length(x) , N+1]);
p0 = zeros([length(x) 1]);
p1 = zeros([length(x) 1]);

p0(:) = ones(xsize)/sqrt(sqrt(pi));
p1(:) = sqrt(2)*x.*p0(:);

p(:,1) = p0;
p(:,2) = p1;

% Loops over the recurrence
for q = 2:N
    p(:,q+1) = sqrt(2)*x.*p(:, q) - sqrt(q-1)*p(:,q-1);
    p(:,q+1) = p(:,q+1)/sqrt(q);
end

return