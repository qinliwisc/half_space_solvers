function [V,D,A_2d,B_2d]=GeneralizedEigen_2d(N,u, dampcoef)

if (nargin == 0) 
    N = 4; u = 0; dampcoef = 0.1;
end


% Calculating A 

[alpha,beta] = half_hermite_shift_recurrence_symb(N-1,0);
Ap = diag(alpha) + diag(sqrt(beta(2:end)), 1) + diag(sqrt(beta(2:end)), -1);
Am = Ap;

A = [Ap(1:end-1,1:end-1)-Am(1:end-1,1:end-1),Ap(1:end-1,:)+Am(1:end-1,:)];
A = [A;Ap(:,1:end-1)+Am(:,1:end-1),Ap-Am];
A = A/2; 

A_2d = kron(eye(size(A)), A);

% The weight is \sqrt{M(u)} \sqrt{M(u+v)}
% Basis is \sqrt{M} polymonials; shifted to centered at -u
[Xp, Wp] = half_hermquad_shift_symb(N-1, u/2);
Wp = exp(-u^2/4) / sqrt(pi) * Wp;
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 

Xp = Xp - u; % shift to center at -u;

[Xm, Wm] = half_hermquad_shift_symb(N-1, -u/2);
Wm = exp(-u^2/4) / sqrt(pi) * Wm;
[Hm] = half_hermite_shift_poly(Xm, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2);

Xm = - Xm - u; % shift to center at -u;

X = [Xp; Xm]; 
W = [Wp; Wm]; 
H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
vu = X+u;

H_2d = kron(H,H); W_2d = kron(W,W);
[yy,xx] = meshgrid(X,X);
vu_2d = kron(ones(size(X)),vu);

% Note that the basis H is not orthonormal with respect to W 
% It is the shifted version of the half space Hermite polynomial 
% Be careful when doing the boundary condition !!!


% Figuring out Xp, Xm, X0 etc.
chi_0x = inline('(vx.^2+vy.^2-2)/sqrt(2)','vx','vy');
chi_0y = inline('sqrt(2)*vy','vx','vy');
chi_p = inline('(vx.^2+vy.^2+2*vx)/2','vx','vy');
chi_m = inline('(vx.^2+vy.^2-2*vx)/2','vx','vy');

chipX = chi_p(xx,yy); chipX = reshape(chipX,length(X)^2,1);
chimX = chi_m(xx,yy); chimX = reshape(chimX,length(X)^2,1);
chi0x = chi_0x(xx,yy); chi0x = reshape(chi0x,length(X)^2,1);
chi0y = chi_0y(xx,yy); chi0y = reshape(chi0y,length(X)^2,1);

if (abs(u)<eps)
    Xp = chipX;
    Xm = chimX;
    X0x = chi0x;
    X0y = chi0y;

    matXp = horzcat(Xp);
    matXm = horzcat(Xm);
    matX0 = horzcat(X0x,X0y);
elseif (u>0) 
    if (abs(u-sqrt(3/2))<eps)
        Xp2 = chipX; 
        Xp1x = chi0x;
        Xp1y = chi0y;
        X0 = chimX;

        matXp = horzcat(Xp2, Xp1x, Xp1y);
        matXm = zeros(size(X0));
        matX0 = horzcat(X0);
    elseif (u<sqrt(3/2))
        Xp2 = chipX; 
        Xp1x = chi0x; 
        Xp1y = chi0y;
        Xm1 = chimX;
        
        matXp = horzcat(Xp2, Xp1x, Xp1y);
        matXm = horzcat(Xm1);
        matX0 = zeros(size(Xm1));
    elseif (u>sqrt(3/2))
        Xp3 = chipX; 
        Xp2x = chi0x; 
        Xp2y = chi0y;
        Xp1 = chimX;
        
        matXp = horzcat(Xp3, Xp2x, Xp2y, Xp1);
        matXm = zeros(size(Xp1));
        matX0 = zeros(size(Xp1));
    end
elseif (u<0)
    if (abs(u+sqrt(3/2))<eps)
        X0 = chipX; 
        Xm2 = chimX;
        Xm1x = chi0x;
        Xm1y = chi0y;
        
        matXp = zeros(size(X0));
        matXm = horzcat(Xm2, Xm1x, Xm1y);
        matX0 = horzcat(X0);
    elseif (u>-sqrt(3/2))
        Xp1 = chipX; 
        Xm1x = chi0x;
        Xm1y = chi0y;
        Xm2 = chimX;
        
        matXp = horzcat(Xp1);
        matXm = horzcat(Xm2, Xm1x, Xm1y);
        matX0 = zeros(size(Xm1x));
    elseif (u<-sqrt(3/2))
        Xm1 = chipX; 
        Xm2x = chi0x;
        Xm2y = chi0y;
        Xm3 = chimX;
        
        matXp = zeros(size(Xm1));
        matXm = horzcat(Xm3, Xm2x, Xm2y, Xm1);
        matX0 = zeros(size(Xm1));
    end
end

% Calculating B

LH_coef_p = matXp.'*diag(W_2d)*H_2d;
LH_coef_m = matXm.'*diag(W_2d)*H_2d;
LH_coef_0 = matX0.'*diag(W_2d)*H_2d;

matvuXp = (vu_2d * ones(1, size(matXp, 2))) .* matXp;
matvuXm = (vu_2d * ones(1, size(matXm, 2))) .* matXm;
matvuX0 = (vu_2d * ones(1, size(matX0, 2))) .* matX0;
matvu_Linv_vuX0 = (vu_2d * ones(1, size(matX0, 2))) .* matvuX0;

LH_coef_vuXp = matvuXp.'*diag(W_2d)*H_2d;
LH_coef_vuXm = matvuXm.'*diag(W_2d)*H_2d;
LH_coef_vuX0 = matvuX0.'*diag(W_2d)*H_2d;
LH_coef_vu_Linv_vuX0 = matvu_Linv_vuX0.'*diag(W_2d)*H_2d;


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


        
B_2d = H_2d.'*diag(W_2d)*LH;
B_2d = (B_2d + B_2d.')/2;
B_2d =  - B_2d - eye(size(B_2d));

[V,D]=eig(A_2d,B_2d);
D = diag(D);

fprintf('Number of negative eigenvalue = %ld (should equal to %ld)\n', ...
        sum(D < -1e-13), (2*N-1)*(N-1));

fprintf('Number of positive eigenvalue = %ld (should equal to %ld)\n', ...
        sum(D > 1e-13), (2*N-1)*(N-1));

fprintf('Number of zero eigenvalue = %ld (should equal to %1d)\n', ...
        sum(abs(D) < 1e-13), 2*N-1);

return