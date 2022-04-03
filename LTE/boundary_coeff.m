function [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = boundary_coeff(N,fb_ini_s)

%% preparation
if (nargin == 0)
    clear;
    N = 80;
    fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
	fb_ini_s.fun = inline('v');
end


%% discretization on v: inner product quadrature
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp/2;

fb_ini_s_data = feval(fb_ini_s.fun,Xp);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

X = [Xp; Xm]; W = [Wp; Wm];

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

%% recovry
chi_data = ones(length(X),1);

[fb_c] = BoundaryData_c(N,0,fb_ini_s); fb = H*fb_c;% plot(X,fb);pause;
[gb_c] = BoundaryData_c(N,1,fb_ini_s); gb = H*gb_c;% plot(X,gb);pause;

C = X'*diag(W)*gb;
D = X'*diag(W)*fb;
eta = D/C;

g_tilde = gb*eta;
phi = chi_data*eta;

Boundary_data = fb - g_tilde + phi;

[val, error] = Hfun(Xp, 1, 0.0001);
real_soln = val/sqrt(3) - Xp;
real_soln_neg = real_soln(end:-1:1);
real_soln = [Xp;real_soln_neg];


figure(1)
plot(X, fb, '.-', X, gb, 'ro-', X, phi, 'gx-');
figure(2)
set(gca,'fontsize',20);
plot(X,Boundary_data,'o',X,real_soln);%Xp,fb_ini_s_data);
% plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
legend('numerical solution','exact solution','Location','NorthWest');
xlabel('v');
ylabel('f_b');
title('f_b (v>0) = v, N = 80');
print(gcf,'-depsc2','v.eps');