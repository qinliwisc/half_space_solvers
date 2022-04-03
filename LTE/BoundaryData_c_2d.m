function [coeff,coeff_2d] = BoundaryData_c(N,index,fb_ini_s,x_space)

% preparation
if (nargin == 0)
    clear;
    N = 36;

    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
    Hp = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp(ind,:);

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm([end:-1:1],:);

    X = [Xp; Xm]; 
    W = [Wp; Wm]; 

    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

    fb_ini_s.fun = inline('v');
    index = 0;
    
    x_space = [0:0.01:5];
end

dampcoef = 0.01;

%% discretization on v
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
Hp = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp/2;

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
X = [Xp;Xm]; W = [Wp;Wm];

%% initial data
fb_ini(:,1) = feval(fb_ini_s.fun,Xp);
fb_ini(:,2) = ones(length(Xp),1);

[V,D,Amat,Bmat] = GeneralizedEigen(N,dampcoef,X,W,H);

if (sum(D>-1e-14) ~= N)
    error('Wrong number of negative eigenvalues!');
end
index_neg = find(D<-1e-14);
Vcons = V(:,find(D>-1e-14));
Vneg = V(:,index_neg);
Dneg = D(index_neg);
% Dneg = zeros(size(D)); Dneg(index_neg) = D(index_neg);

% consistency check for the eignevectors
if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
    error('Something wrong in the generalized eigenvalue');
end

A = zeros(2*N-1, 2*N-1);
b = zeros(2*N-1, 1);

% The first N rows are constraints from the growing and zero modes
A(1:N, :) = Vcons';
b(1:N, :) = 0;

% The rest rows are the Galerkin condition from the boundary values
A(N+1:end, :) = Hp(:, 1:N-1)' * diag(Wp.*Xp) * [Hp(:, 1:N-1), Hp];
b(N+1:end, :) = Hp(:, 1:N-1)' * diag(Wp.*Xp) * fb_ini(:,index+1);
% A(N+1:end, :) = Hp(:, 1:N-1).' * diag(Wp.*Xp) * [Hp(:, 1:N-1), Hp];
% b(N+1:end, :) = Hp(1:N, 1:N-1).' * diag(Wp.*Xp) * fb_ini(:,index+1);

coeff = A \ b;

thresholding = zeros(2*N-1,1);
thresholding(1:N-1) = [1:1:N-1];
thresholding = thresholding/N;
thresholding = (1+cos(pi*thresholding))/2;

index = ones(length(coeff),1);
index = (abs(coeff)>2.5*1e-3);

coeff = coeff.*index;

coeff_x = (V'*coeff);
coeff_x_2d = coeff_x * ones(1,length(x_space));
spatial = ones(2*N-1,length(x_space));
spatial(index_neg,:) = exp(1./Dneg*x_space);
coeff_x_2d = coeff_x_2d.*spatial;
coeff_2d = (V')\coeff_x_2d;

% coeff_x = (Vneg'*coeff);
% coeff_x_2d = coeff_x * ones(1,length(x_space));
% spatial = exp(1./Dneg*x_space);
% coeff_x_2d = coeff_x_2d.*spatial;
% coeff_2d = Vneg*coeff_x_2d;

return