function coeff = BoundaryData_multi(N,index,fb_ini_s1,fb_ini_s2,sigma)

% preparation
if (nargin == 0)
    clear;
    N = 80;

    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
    Hp = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp(ind,:);

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm([end:-1:1],:);

    X = [Xp; Xm]; 
    W = [Wp; Wm]; 

    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

    fb_ini_s1.fun = inline('v');
    fb_ini_s2.fun = inline('2*v');
    
    index = 0;
    
    sigma11 = 0.5; sigma12 = 0.5; sigma21 = 0.5; sigma22 = 0.5;
%     sigma11 = 1; sigma12 = 0; sigma21 = 0; sigma22 = 1;
%     sigma11 = 0.25; sigma12 = 0.25; sigma21 = 0.75; sigma22 = 0.75;
    sigma = [sigma11,sigma12;sigma21,sigma22];
end

dampcoef = 0.01;

%% discretization on v
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
Hp = half_legendre_poly(Xp,[0:1:N-1]);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
X = [Xp;Xm]; W = [Wp;Wm]; W = W/2;

%% initial data
fb_ini(:,1) = [feval(fb_ini_s1.fun,Xp);feval(fb_ini_s2.fun,Xp)];
fb_ini(:,2) = [sigma(1,2)*ones(length(Xp),1);(1-sigma(1,1))*ones(length(Xp),1)];

[V,D,Amat,Bmat] = GeneralizedEigen_multi(N,dampcoef,X,W,H,sigma);

if (sum(D>-1e-14) ~= 2*N)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

% consistency check for the eignevectors
% if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
%     error('Something wrong in the generalized eigenvalue');
% end

A = zeros(4*N-2, 4*N-2);
b = zeros(4*N-2, 1);

% The first N rows are constraints from the growing and zero modes
A(1:2*N, :) = Vcons.'*Bmat; 
b(1:2*N, :) = 0;

% The rest rows are the Galerkin condition from the boundary values
% A(2*N+1:end, :) = Hp(:, 1:N-1)' * diag(Wp.*Xp) * [Hp(:, 1:N-1), Hp];
% b(2*N+1:end, :) = Hp(:, 1:N-1)' * diag(Wp.*Xp) * fb_ini(:,index+1);
% A(N+1:end, :) = Hp(:, 1:N-1).' * diag(Wp.*Xp) * [Hp(:, 1:N-1), Hp];
% b(N+1:end, :) = Hp(1:N, 1:N-1).' * diag(Wp.*Xp) * fb_ini(:,index+1);
term = Hp(1:N, 1:N-1).' * diag(Wp.*Xp);
term = [term,zeros(size(term));zeros(size(term)),term];
term2 = [Hp(:, 1:N-1), Hp];
term2 = [term2,zeros(size(term2));zeros(size(term2)),term2];
A(2*N+1:end, :) = term * term2;
b(2*N+1:end, :) = term * fb_ini(:,index+1);

coeff = A \ b;

thresholding = zeros(2*N-1,1);
thresholding(1:N-1) = [1:1:N-1];
thresholding = thresholding/N;
thresholding = (1+cos(pi*thresholding))/2;

% %index = ones(length(coeff),1);
% index = (abs(coeff)>2.5*1e-3);

coeff(1:2*N-1) = coeff(1:2*N-1).*thresholding;%index;
coeff(2*N:4*N-2) = coeff(2*N:4*N-2).*thresholding;%index;

fb1 = H*(coeff(1:2*N-1));
fb2 = H*(coeff(2*N:4*N-2));
plot(X,fb1,'.-.',X,fb2);
return