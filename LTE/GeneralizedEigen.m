function [V,D,A,B]=GeneralizedEigen(N,dampcoef,X,W,H)

if nargin == 0
    N = 16; dampcoef = 0.01;
    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend');
    [Hp] = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp/2;

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
%     [Hm] = half_legendre_poly(Xm, [0:1:N-1]);
    Hm = Hp; Hm = Hm(ind,:);

    X = [Xp; Xm]; W = [Wp; Wm]; 
    
    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
end

% Calculating A

[alpha,beta] = half_legendre_recurrence(N-1);
Ap = diag(alpha) + diag(beta(1:end-1), 1) + diag(beta(1:end-1), -1);

A = [zeros(N-1,N-1),Ap(1:end-1,:)];
A = [A;Ap(:,1:end-1),zeros(N,N)]/2;

chi0X = ones(length(X),1);
matX0 = chi0X(:);

vu = X;

% Calculating B
LH_coef_0 = matX0.'*diag(W)*H;

matvuX0 = (vu * ones(1, size(matX0, 2))) .* matX0;
matvu_Linv_vuX0 = (vu * ones(1, size(matX0, 2))) .* matvuX0;

LH_coef_vuX0 = matvuX0.'*diag(W)*H;
LH_coef_vu_Linv_vuX0 = matvu_Linv_vuX0.'*diag(W)*H;


% BGK collision
LH = - matX0 * LH_coef_0; 
% damp X_0
LH = LH + dampcoef * matvuX0 * LH_coef_vuX0;
% damp L^{-1}((v+u) X_0)
LH = LH + dampcoef * matvu_Linv_vuX0 * LH_coef_vu_Linv_vuX0;

B = H.'*diag(W)*LH;
B = (B + B.')/2;
B =  - B - eye(size(B));

[V,D]=eig(A,B);
D = diag(D);

fprintf('Number of negative eigenvalue = %ld (should equal to %ld)\n', ...
        sum(D < -eps), N-1);

% number of negative eigenvalues equals to number of even basis
% functions

fprintf('Number of zero eigenvalue = %ld (should equal to 1)\n', ...
        sum(abs(D) < 1e-13));

return