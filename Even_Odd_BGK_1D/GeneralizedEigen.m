function [V,D,A,B]=GeneralizedEigen(N,u, dampcoef,X,W,H)

vu = X+u;

% Calculating A 

[alpha,beta] = half_hermite_shift_recurrence_symb(N-1,0);
Ap = diag(alpha) + diag(sqrt(beta(2:end)), 1) + diag(sqrt(beta(2:end)), -1);
Am = Ap; 

A = [Ap(1:end-1,1:end-1)-Am(1:end-1,1:end-1),Ap(1:end-1,:)+Am(1:end-1,:)];
A = [A;Ap(:,1:end-1)+Am(:,1:end-1),Ap-Am];
A = A/2; 

chi_0 = inline('(2*v.^2-3)/sqrt(6)');
chi_p = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
chi_m = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');

chipX = chi_p(X); 
chimX = chi_m(X); 
chi0X = chi_0(X); 

if (abs(u)<eps)
    Xp = chipX; 
    Xm = chimX;
    X0 = chi0X;

    matXp = horzcat(Xp);
    matXm = horzcat(Xm);
    matX0 = horzcat(X0);
elseif (u>0) 
    if (abs(u-sqrt(3/2))<eps)
        Xp2 = chipX; 
        Xp1 = chi0X;
        X0 = chimX;

        matXp = horzcat(Xp2, Xp1);
        matXm = zeros(size(X0));
        matX0 = horzcat(X0);
    elseif (u<sqrt(3/2))
        Xp2 = chipX; 
        Xp1 = chi0X;
        Xm1 = chimX;
        
        matXp = horzcat(Xp2, Xp1);
        matXm = horzcat(Xm1);
        matX0 = zeros(size(Xm1));
    elseif (u>sqrt(3/2))
        Xp3 = chipX; 
        Xp2 = chi0X;
        Xp1 = chimX;
        
        matXp = horzcat(Xp3, Xp2, Xp1);
        matXm = zeros(size(Xp1));
        matX0 = zeros(size(Xp1));
    end
elseif (u<0)
    if (abs(u+sqrt(3/2))<eps)
        X0 = chipX; 
        Xm2 = chimX;
        Xm1 = chi0X;
        
        matXp = zeros(size(X0));
        matXm = horzcat(Xm2, Xm1);
        matX0 = horzcat(X0);
    elseif (u>-sqrt(3/2))
        Xp1 = chipX; 
        Xm1 = chi0X;
        Xm2 = chimX;
        
        matXp = horzcat(Xp1);
        matXm = horzcat(Xm2, Xm1);
        matX0 = zeros(size(Xm1));
    elseif (u<-sqrt(3/2))
        Xm1 = chipX; 
        Xm2 = chi0X;
        Xm3 = chimX;
        
        matXp = zeros(size(Xm1));
        matXm = horzcat(Xm3, Xm2, Xm1);
        matX0 = zeros(size(Xm1));
    end
end

% Calculating B

LH_coef_p = matXp.'*diag(W)*H;
LH_coef_m = matXm.'*diag(W)*H;
LH_coef_0 = matX0.'*diag(W)*H;

matvuXp = (vu * ones(1, size(matXp, 2))) .* matXp;
matvuXm = (vu * ones(1, size(matXm, 2))) .* matXm;
matvuX0 = (vu * ones(1, size(matX0, 2))) .* matX0;
matvu_Linv_vuX0 = (vu * ones(1, size(matX0, 2))) .* matvuX0;

LH_coef_vuXp = matvuXp.'*diag(W)*H;
LH_coef_vuXm = matvuXm.'*diag(W)*H;
LH_coef_vuX0 = matvuX0.'*diag(W)*H;
LH_coef_vu_Linv_vuX0 = matvu_Linv_vuX0.'*diag(W)*H;


% BGK collision
LH = - matXp * LH_coef_p - matXm * LH_coef_m - matX0 * LH_coef_0; 
% damp X_+ 
LH = LH + dampcoef * matvuXp * LH_coef_vuXp;
% damp X_-
LH = LH + dampcoef * matvuXm * LH_coef_vuXm;

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