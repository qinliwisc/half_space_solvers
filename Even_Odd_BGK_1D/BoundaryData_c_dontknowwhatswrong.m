function [coeff]=BoundaryData_c(N,u,index,X_q,W_q,H_q,fb_ini_s)

if (nargin ==0)
    N = 16; u = 0.5; index = 1;
    [Xp, Wp] = half_hermquad_shift_symb(N-1, u/2);
    Wp = exp(-u^2/4) / sqrt(pi) * Wp;
    [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
    [Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 

    Xp = Xp - u; % shift to center at -u;

    [Xm, Wm] = half_hermquad_shift_symb(N-1, -u/2);
    Wm = exp(-u^2/4) / sqrt(pi) * Wm;
    [Xm,ind] = sort(Xm,'ascend'); Wm = Wm(ind);
    [Hm] = half_hermite_shift_poly(Xm, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2);

    Xm = - Xm - u; % shift to center at -u;

    X_q = [Xp; Xm]; 
    W_q = [Wp; Wm]; 

    H_q = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
    fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');
end

%% preparation
dampcoef = 0.1;

coeff = [];

chi_0 = inline('(2*v.^2-3)/sqrt(6)');
chi_p = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
chi_m = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(N-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp-u;
Xm = Xm-u; 


X = [Xp;Xm]; W = [Wp;Wm];
vu = X + u;

%% initial data
% Note that the initial data is fb_ini \sqrt{M(v+u)} !!
fb_ini(:,1) = feval(fb_ini_s.fun,Xp + u);
fb_ini(:,2) = chi_0(Xp + u);
fb_ini(:,3) = chi_p(Xp + u);
fb_ini(:,4) = chi_m(Xp + u);

%% 
H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

[V,D,Amat,Bmat] = GeneralizedEigen(N,u,dampcoef,X_q,W_q,H_q);
if (sum(D>-1e-14) ~= N)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

% consistency check for the eignevectors
if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
    error('Something wrong in the generalized eigenvalue');
end

%% linear system for c
A = zeros(2*N-1,2*N-1);
b = zeros(2*N-1, 1);

% The first N rows are constraints from the growing and zero modes
A(1:N, :) = Vcons.'*Bmat; 
b(1:N, :) = 0;

for round = 1:length(index)
    % The rest rows are the Galerkin condition from the boundary values
    A(N+1:end, :) = Hp(:, 1:N-1).' * diag(Wp.*(Xp + u)) * [Hp(:, 1:N-1), Hp];
    b(N+1:end, :) = Hp(:, 1:N-1).' * diag(Wp.*(Xp + u)) * fb_ini(:,index(round)+1);

    c = A \ b; 
    coeff = horzcat(coeff,c);

    % consistency check for the solution 
    if (norm(b - A * c)>1e-8)
        error('Something wrong in the linear solver');
    end
end

return