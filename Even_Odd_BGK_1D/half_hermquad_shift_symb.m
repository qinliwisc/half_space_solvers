function [X,W]=half_hermquad_shift_symb(n,u) % the method presented in Golub's paper
if nargin==0
    n = 15; u = 0;
end

[alpha,beta] = half_hermite_shift_recurrence_symb(n,u);

pi_s = sym('pi');
if (abs(u) < eps)
    u_s = 0;
elseif (abs(u - sqrt(3/2)) < eps)
    u_s = sym(sqrt(3/2));
elseif (abs(u + sqrt(3/2)) < eps)
    u_s = sym(-sqrt(3/2));
else
    u_s = sym(u);
end
temp = 1 + erf(u_s);
temp_s = sym(temp);
temp_s = temp_s * sqrt(pi_s)/2;
% temp2_s = exp(-u_s^2)/2;
% 
% alpha = sym(zeros(n, 1));
% beta = sym(zeros(n, 1));
% 
% 
% alpha(1) = temp2_s + u_s*temp_s;
% %alpha(1) = temp2_s + u*temp_s;
% alpha(1) = alpha(1)/temp_s;
% % alpha(2) = 2/(pi_s-2)/sqrt(pi_s);
% % beta(2) = (pi_s-2)/pi_s/2;%alpha(1)/(alpha(1)+alpha(2))/2;
% 
% alpha = vpa(alpha);
% beta = vpa(beta);
% 
% for k=1:n-1
%     beta(k+1) = (k-1) + 1/2 - alpha(k)^2-beta(k);
%     beta(k+1) = beta(k+1) + u_s*alpha(k);
%     alpha(k+1) = sum(alpha(1:k))/2/beta(k+1) - alpha(k);
%     alpha(k+1) = alpha(k+1) + u_s;
% end


A = diag(alpha(1:end)) + diag(sqrt(beta(2:end)), 1) + ...
    diag(sqrt(beta(2:end)), -1);
%A = (A + A.')/2;
[V,X]=eig(A);
X = diag(X);
W = vpa(V(1, :)'.^2) * temp_s;

X = double(X);
W = double(W);