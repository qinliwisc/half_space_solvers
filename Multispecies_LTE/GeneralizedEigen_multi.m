function [V,D,A,B]=GeneralizedEigen_multi(N,dampcoef,X,W,H,sigma)

if nargin == 0
    N = 8; dampcoef = 0.01;
    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend');
    [Hp] = half_legendre_poly(Xp,[0:1:N-1]);

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm(ind,:);

    X = [Xp; Xm]; W = [Wp; Wm]; W = W/2;
    
    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
    
    sigma11 = 0.5; sigma12 = 0.5; sigma21 = 0.5; sigma22 = 0.5;
    sigma = [sigma11,sigma12;sigma21,sigma22];
end

H_tilde = [H,zeros(size(H));zeros(size(H)),H];
Hv = (X*ones(1,2*N-1)).*H;
Hv_tilde = [Hv,zeros(size(Hv));zeros(size(Hv)),Hv];

% Calculating A
A = (H_tilde)'*diag([W;W])*Hv_tilde;
A(abs(A)<1e-14) = 0;

% Calculating B
term1 = zeros(2*(2*N-1),1); term1(1) = 1;
B = term1*[sigma(1,1)*W',sigma(1,2)*W']*H_tilde;
term2 = zeros(2*(2*N-1),1); term2(2*N) = 1;
B = B + term2*[sigma(2,1)*W',sigma(2,2)*W']*H_tilde;
B = B - eye(2*(2*N-1));
term3 = sigma(1,2)*[H'*diag(W)*X;zeros(2*N-1,1)];
term4 = (1-sigma(1,1))*[zeros(2*N-1,1);H'*diag(W)*X];
alpha = [sigma(1,2)*X',(1-sigma(1,1))*X']*diag([W;W])*H_tilde;
B = B - dampcoef*(term3+term4)*alpha;
term5 = sigma(1,2)*[H'*diag(W)*(X.^2);zeros(2*N-1,1)];
term6 = (1-sigma(1,1))*[zeros(2*N-1,1);H'*diag(W)*(X.^2)];
beta = [sigma(1,2)*(X.^2)',(1-sigma(1,1))*(X.^2)']*diag([W;W])*H_tilde;
B = B - dampcoef*(term5+term6)*beta;

[V,D]=eig(A,B);
D = diag(D);

fprintf('Number of negative eigenvalue = %ld (should equal to %ld)\n', ...
        sum(D < -eps), 2*(N-1));

% number of negative eigenvalues equals to number of even basis
% functions

fprintf('Number of zero eigenvalue = %ld (should equal to 2)\n', ...
        sum(abs(D) < 1e-13));

return