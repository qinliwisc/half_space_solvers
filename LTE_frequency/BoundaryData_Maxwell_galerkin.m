function coeff = BoundaryData_Maxwell_galerkin(N,N_omega,index,alpha1,alpha2,alpha3,fb_ini)

% preparation
if (nargin == 0)
    clear;
    N_omega = 8; N = 16; index = 1;
    fb_ini.fun = inline('v - v + omega.^(-1/2)','v','omega');
    alpha1 = 0.3; alpha2 = 0.3; alpha3 = 0.4;
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

H_even = kron(diag(sqrt(beta_inv)),H_X(:,1:N-1));
H_odd = kron(diag(sqrt(beta_inv)),H_X(:,N:end));

H_even_p = kron(diag(sqrt(beta_inv)),Hp(:,1:end-1));
Weight_pp = kron(beta,Wp/2); Weight_pp = Weight_pp / sum(Weight);
Weight_mm = kron(beta,Wm/2); Weight_mm = Weight_mm / sum(Weight);
X_pp = kron(ones(N_omega,1),Xp);
X_mm = kron(ones(N_omega,1),Xm);

Equi_p = kron(alpha_inv,ones(length(Xp),1));
Equi_m = kron(alpha_inv,ones(length(Xm),1));

% Hpp = kron(diag(sqrt(beta_inv)),Hp(:,1:end-1));
% Hpp = Hpp/sqrt(iii(1,1));
% Weight_pp = kron(beta,Wp/2); Weight_pp = Weight_pp / sum(Weight);
% Weight_mm = kron(beta,Wm/2); Weight_mm = Weight_mm / sum(Weight);
% H_pos = kron(diag(sqrt(beta_inv)),[Hp(:,1:end-1),Hp]);
% H_pos = H_pos/sqrt(iii(1,1));
% X_pp = kron(ones(N_omega,1),Xp);
% X_mm = kron(ones(N_omega,1),Xm);
% 
% 
% H_neg_updown = kron(diag(sqrt(beta_inv)),[Hm_updown(:,1:end-1),Hm_updown]);
% H_neg = kron(diag(sqrt(beta_inv)),[Hm(:,1:end-1),Hm]);
% H_neg_updown = H_neg_updown/sqrt(iii(1,1));
% H_neg = H_neg/sqrt(iii(1,1));
% 
% Equi_p = kron(alpha_inv,ones(length(Xp),1));
% Equi_m = kron(alpha_inv,ones(length(Xm),1));

%% discretization on v
[omegaomega,xx] = meshgrid(omega,X);

%% initial data
if index == 0
    fb_ini = feval(fb_ini.fun,xx(1:length(Xp),:),omegaomega(1:length(Xp),:));
    fb_ini_1d = reshape(fb_ini,N_omega*length(Xp),1);
elseif index == 1
    fb_ini1 = Equi_p;
    fb_ini1 = reshape(fb_ini1,N_omega*length(Xp),1);
    fb_ini2 = Equi_m;
    fb_ini2 = reshape(fb_ini2,N_omega*length(Xm),1);
    fb_ini3 = Equi_m;
    fb_ini3 = reshape(fb_ini3,N_omega*length(Xm),1);
    fb_ini3 = Equi_p * (X_mm.*Equi_m)'*diag(Weight_mm)*fb_ini3;
    fb_ini4 = (X_pp)'*diag(Weight_pp)*(Equi_p.^2);
    fb_ini3 = fb_ini3/fb_ini4;
    fb_ini_1d = fb_ini1 - alpha2 * fb_ini2 + alpha3*fb_ini3;
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
% term1 = Hpp' * diag(Weight_pp.*X_pp) * H_pos;
% term2 = Hpp' * diag(Weight_pp.*X_pp) * H_neg_updown;
% term3 = Hpp' * diag(Weight_pp.*X_pp) * Equi_p * (X_mm.*Equi_m)'*diag(Weight_mm)*H_neg;
% term4 = (X_pp)'*diag(Weight_pp)*(Equi_p.^2);
% term3 = term3/term4;
% term = term1 - alpha2 * term2 + alpha3 * term3;
% A(N*N_omega+1:end, :) = term;
% b(N*N_omega+1:end, :) = alpha1 * Hpp' * diag(Weight_pp.*X_pp) * fb_ini_1d;
term = H_even' * diag(Weight.*X_2d) * H_odd;
A(N*N_omega+1:end,(N-1)*N_omega+1:end) = H_even' * diag(Weight.*X_2d) * H_odd;
term = Equi_p' * diag(Weight_pp.*X_pp) * Equi_p;
term = 1/term*Equi_p;
term2 = Equi_p' * diag(Weight_pp.*X_pp) * H_even_p;
term2 = term*term2;
term3 = H_even_p - 2*alpha2/(1-alpha3)/(1+alpha2+alpha3)*term2;
term3 = (1-alpha3)/(1+alpha3)*term3;
term = 2*H_even_p' * diag(Weight_pp.*X_pp) * term3;
A(N*N_omega+1:end,1:(N-1)*N_omega) = term;

% A(N*N_omega+1:end,1:(N-1)*N_omega) = 2*Equi_p' * diag(Weight_pp.*X_pp) * term3;

coeff = A \ b;

f_recover = H*coeff;
f_recover = reshape(f_recover,length(X),N_omega);
return