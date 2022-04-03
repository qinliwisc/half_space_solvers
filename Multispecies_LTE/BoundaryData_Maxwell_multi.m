function coeff = BoundaryData_Maxwell_multi(N,index,fb_ini_s1,fb_ini_s2,alpha,sigma)

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
    fb_ini_s2.fun = inline('v');
    
    index = 0;
    
    alpha = [0.1,0.9/2,0.9/2];
%     alpha = [1,0,0];
    sigma11 = 0.5; sigma12 = 0.5; sigma21 = 0.5; sigma22 = 0.5;
%     sigma11 = 1; sigma12 = 0; sigma21 = 0; sigma22 = 1;
%     sigma11 = 0.25; sigma12 = 0.25; sigma21 = 0.75; sigma22 = 0.75;
    sigma = [sigma11,sigma12;sigma21,sigma22];
end

dampcoef = 0.01;

%% discretization on v
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind); Wp = Wp/2;
Hp = half_legendre_poly(Xp,[0:1:N-1]);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
X = [Xp;Xm]; W = [Wp;Wm];

%% initial data
fb_ini(:,1) = [feval(fb_ini_s1.fun,Xp);feval(fb_ini_s2.fun,Xp)];
fb_ini(:,2) = [sigma(1,2)*ones(length(Xp),1);(1-sigma(1,1))*ones(length(Xp),1)];

%% negative mode preparation
[V,D,Amat,Bmat] = GeneralizedEigen_multi(N,dampcoef,X,W,H,sigma);

if (sum(D>-1e-14) ~= 2*N)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

%% boundary condition
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

% term = H(1:length(Xp),:) - 2*alpha(2)*ones(length(Xp),1)*Xm'*diag(Wm)*H(length(Xp)+1:end,:) - alpha(3)*flipud(H(1:length(Xp),:));
% term = Hp(:,1:end-1)'*diag(Wp)*term;
% term = [term,zeros(size(term));zeros(size(term)),term];
H_pos = H(1:length(Xp),:); H_neg = H(length(Xp)+1:end,:);
term1 = [H_pos,zeros(size(H_pos));zeros(size(H_pos)),H_pos];
term3 = alpha(3)*[flipud(H_neg),zeros(size(H_neg));zeros(size(H_neg)),flipud(H_neg)];
term2 = fb_ini(:,2)*[Xm;Xm]'*diag([Wm;Wm])*[H_neg,zeros(size(H_neg));zeros(size(H_neg)),H_neg];
term2 = -2*alpha(2)*term2/(sigma(1,2)+1-sigma(1,1));
term = term1 - term2 - term3;

term4 = [Hp(:,1:end-1),zeros(size(Hp(:,1:end-1)));zeros(size(Hp(:,1:end-1))),Hp(:,1:end-1)];

A(2*N+1:end, :) = term4'*diag([Wp;Wp])*term;
b(2*N+1:end, :) = alpha(1)*term4'*diag([Wp;Wp])*fb_ini(:,index+1);

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