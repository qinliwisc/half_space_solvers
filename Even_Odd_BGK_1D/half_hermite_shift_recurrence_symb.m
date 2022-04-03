function [alpha,beta] = half_hermite_shift_recurrence_symb(N,u)
% recurrence -- recurrence coefficients for Hermite polynomials
%
% [a,b] = recurrence(N)
%     Calculates the first N recurrence coefficients for the hermite
%     polynomials.  


if nargin == 0
    N = 15; u = sqrt(3/2);
end

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
m_0_s = 1 + erf(u_s);
m_0_s = m_0_s * sqrt(pi_s)/2;
m_1_s = exp(-u_s^2)/2 + u_s*m_0_s;
m_2_s = -u_s*exp(-u_s^2)/2 + m_0_s/2 + 2*u_s*m_1_s - u_s^2*m_0_s;

alpha = sym(zeros(N+1, 1));
beta = sym(zeros(N+1, 1));


alpha(1) = m_1_s/m_0_s;
beta(2) = m_2_s-2*alpha(1)*m_1_s+alpha(1)^2*m_0_s;
beta(2) = beta(2)/m_0_s;
alpha(2) = alpha(1)/2/beta(2) - alpha(1);
alpha(2) = alpha(2) + u_s;

alpha = vpa(alpha, 64);
beta = vpa(beta, 64);

for k=2:N
    beta(k+1) = (k-1) + 1/2 - alpha(k)^2-beta(k);
    beta(k+1) = beta(k+1) + u_s*alpha(k);
    alpha(k+1) = sum(alpha(1:k))/2/beta(k+1) - alpha(k);
    alpha(k+1) = alpha(k+1) + u_s;
    alpha = vpa(alpha, 64);
    beta = vpa(beta, 64);
end

alpha = double(alpha);
beta = double(beta);
return