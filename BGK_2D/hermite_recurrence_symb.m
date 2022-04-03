function [alpha] = hermite_recurrence_symb(N)

if nargin == 0
    N = 15;
end

alpha = sym(zeros(N+1,1));
for k = 1:N
    alpha(k+1) = sym(sqrt(k));
end

alpha = double(alpha)/sqrt(2);
return