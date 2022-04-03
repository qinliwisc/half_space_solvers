function [V,D,A,B]=GeneralizedEigen(N,N_omega,dampcoef,X,X_2d,alpha,beta,Weight,H)

if nargin == 0
    N = 8; N_omega = 8; dampcoef = 0.01;
    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend');
    [Hp] = half_legendre_poly(Xp,[0:1:N-1]);

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm(ind,:);

    H_X = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
    X = [Xp;Xm]; W = [Wp;Wm]; W = W/2;

    omega = [1:1:N_omega];
    alpha = ones(length(omega),1); alpha = alpha(:);
    beta = omega.*exp(-omega/1000); beta = beta(:);

    beta_inv = ones(size(beta))./beta;
    H = kron(diag(sqrt(beta_inv)),H_X);
    Weight = kron(beta,W); Weight = Weight / sum(Weight);
    iii = H'*diag(Weight)*H;
    H = H/sqrt(iii(1,1));
    X_2d = kron(ones(N_omega,1),X);
end

% Calculating A
A = H'*diag(X_2d.*Weight)*H;
A(abs(A)<1e-14) = 0;

% Calculating B
alpha_stretch = alpha*ones(1,length(X));
alpha_stretch = reshape(alpha_stretch',N_omega*length(X),1);
beta_stretch = beta*ones(1,length(X));
beta_stretch = reshape(beta_stretch',N_omega*length(X),1);

term1 = alpha_stretch'*diag(Weight)*H;
term1 = term1'*term1;

term2 = H'*diag(Weight.*alpha_stretch.^2)*H;

term3 = (X_2d./alpha_stretch)'*diag(Weight)*H;
term3 = term3'*term3;
term3 = dampcoef*sum(Weight./(alpha_stretch.^2))*term3;

term4 = (X_2d.^2/alpha_stretch.^3)'*diag(Weight)*H;
term4 = term4'*term4;
term4 = dampcoef*sum(Weight./(alpha_stretch.^2))*term4;

B = term1 - term2 - term3 - term4;
B(abs(B)<1e-14) = 0;

[V,D]=eig(A,B);
D = diag(D);

fprintf('Number of negative eigenvalue = %ld (should equal to %ld)\n', ...
        sum(D < -eps), (N-1)*N_omega);

% number of negative eigenvalues equals to number of even basis
% functions

fprintf('Number of zero eigenvalue = %ld (should equal to %ld)\n', ...
        sum(abs(D) < 1e-13), N_omega);

return