function [c]=BoundaryData_2d(index)

if (nargin == 0)
    N = 8; u = 0; dampcoef = 0.1; index = 0;
end

% 2D weight = exp(-(vx^2+vy^2)/2)
chi_0x = inline('(vx.^2+vy.^2-2)/sqrt(2)','vx','vy');
chi_0y = inline('sqrt(2)*vy','vx','vy');
chi_p = inline('(vx.^2+vy.^2+2*vx)/2','vx','vy');
chi_m = inline('(vx.^2+vy.^2-2*vx)/2','vx','vy');

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(N-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp-u;
Xm = Xm-u; 


X = [Xp;Xm]; W = [Wp;Wm];
vu = X + u;

[yy,xx] = meshgrid(X,vu);
%% initial data
if index == 0
    fb_ini = chi_0x(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(X),1);
elseif index == 1
    fb_ini = chi_0y(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(X),1);
elseif index == 2
    fb_ini = chi_p(xx(:,1:length(Xp)),yy(:,1:length(Xp)));
    fb_ini = reshape(fb_ini,length(Xp)*length(X),1);
elseif index == 3
    fb_ini = chi_m(xx(:,1:length(Xp)),yy(:,1:length(Xp)));
    fb_ini = reshape(fb_ini,length(Xp)*length(X),1);
end

%% polynomials basis
H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
H_2d = kron(H,H);
Wp_2d = kron(W,Wp);
x_2d = kron(ones(size(X)),vu(1:length(Xp)));

Hpp = kron(H,Hp(:,1:end-1)); H_pos = kron(H,[Hp(:,1:end-1),Hp]);

%% linear system
% generalized eigenvalue
[V,D,Amat,Bmat] = GeneralizedEigen_2d(N,u,dampcoef);
if (sum(D>-1e-13) ~= N*(2*N-1))
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

% consistency check for the eignevectors
if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
    error('Something wrong in the generalized eigenvalue');
end

% full system
A = zeros((2*N-1)^2,(2*N-1)^2);
b = zeros((2*N-1)^2, 1);

% The first N(2N-1) rows are constraints from the growing and zero modes
A(1:N*(2*N-1), :) = Vcons.'*Bmat; 
b(1:N*(2*N-1), :) = 0;

%% Collocation
A(N*(2*N-1)+1:end,:) = Hpp'*diag(Wp_2d.*x_2d)*H_pos;
b(N*(2*N-1)+1:end, :) = Hpp'*diag(Wp_2d.*x_2d)*fb_ini;

c = A \ b; 

exponential = exp(-xx.^2/2 - yy.^2/2);

fb = H_2d * c;
fb = reshape(fb,length(X),length(X));
% Fb = fb.*exponential;

Fb_ini = reshape(fb_ini,length(Xp),length(X));
% Fb_ini = Fb_ini.*exponential(1:length(Xp),:);

figure(1)
mesh(Fb)
figure(2)
mesh(Fb_ini - Fb(1:length(Xp),:))
legend('real: u mode','recover','location','southeast');
xlabel('v');
ylabel('f_b');
title('given f_b (v>0) = u mode');

return
