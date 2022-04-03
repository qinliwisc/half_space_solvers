function coeff = BoundaryData(N,N_omega,index,fb_ini)

% preparation
if (nargin == 0)
    clear;
    N_omega = 8; N = 16; index = 0;
    fb_ini.fun = inline('v - v + omega.^(-1/2)','v','omega');
%     fb_ini.fun = inline('v+exp(-omega)','v','omega');
end
dampcoef = 0.01;

[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
Hp = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp(ind,:);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

X = [Xp; Xm]; W = [Wp; Wm]; W = W/2;

H_X = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

omega = [1:1:N_omega];
alpha = omega.^(1/2); alpha = alpha(:);
beta = omega.*exp(-omega/1000); beta = beta(:);

beta_inv = ones(size(beta))./beta;
alpha_inv = ones(size(alpha))./alpha;
H = kron(diag(sqrt(beta_inv)),H_X);
Weight = kron(beta,W); Weight = Weight / sum(Weight);
iii = H'*diag(Weight)*H;
H = H/sqrt(iii(1,1));
X_2d = kron(ones(N_omega,1),X);

Hpp = kron(diag(sqrt(beta_inv)),Hp(:,1:end-1));
Hpp = Hpp/sqrt(iii(1,1));
Weight_pp = kron(beta,Wp/2); Weight_pp = Weight_pp / sum(Weight);
H_pos = kron(diag(sqrt(beta_inv)),[Hp(:,1:end-1),Hp]);
H_pos = H_pos/sqrt(iii(1,1));
X_pp = kron(ones(N_omega,1),Xp);

%% discretization on v
[omegaomega,xx] = meshgrid(omega,X);

%% initial data
if index == 0
    fb_ini = feval(fb_ini.fun,xx(1:length(Xp),:),omegaomega(1:length(Xp),:));
    mesh(omega,Xp,fb_ini);
    fb_ini_1d = reshape(fb_ini,N_omega*length(Xp),1);
end
if index == 1
    fb_ini = ones(length(Xp),1)*alpha_inv';
    mesh(omega,Xp,fb_ini);
    fb_ini_1d = reshape(fb_ini,N_omega*length(Xp),1);
end

[V,D,Amat,Bmat] = GeneralizedEigen(N,N_omega,dampcoef,X,X_2d,...
    alpha,beta,Weight,H);

if (sum(D>-1e-14) ~= N_omega*N)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

A = zeros((2*N-1)*N_omega, (2*N-1)*N_omega);
b = zeros((2*N-1)*N_omega, 1);

% The first N rows are constraints from the growing and zero modes
A(1:N*N_omega, :) = Vcons.'*Bmat; 
b(1:N*N_omega, :) = 0;

% The rest rows are the Galerkin condition from the boundary values
term = Hpp' * diag(Weight_pp.*X_pp);
A(N*N_omega+1:end, :) = term * H_pos;
b(N*N_omega+1:end, :) = term * fb_ini_1d;

coeff = A \ b;

f_recover = H*coeff;
f_recover = reshape(f_recover,length(X),N_omega);

% figure(1);mesh(omega,Xp,fb_ini);
figure(2);mesh(omega,X,f_recover);
% figure(3);mesh(omega,Xp,f_recover(1:length(Xp),:)-fb_ini);
return