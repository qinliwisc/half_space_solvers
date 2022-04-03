function [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = boundary_coeff(N,fb_ini_s)

%% preparation
if (nargin == 0)
    clear;
    N = 80;
    fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
	fb_ini_s.fun = inline('v');
    
    x_space = [0:0.01:0.5]; Nx = length(x_space);
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

[fb_c,fb_c_2d] = BoundaryData_c(N,0,fb_ini_s,x_space); fb = H*fb_c;% plot(X,fb);pause;
[gb_c,gb_c_2d] = BoundaryData_c(N,1,fb_ini_s,x_space); gb = H*gb_c;% plot(X,gb);pause;

C = X'*diag(W)*gb;
D = X'*diag(W)*fb;
eta = D/C;

g_tilde = gb*eta; phi = chi_data*eta;

f_2d = H*fb_c_2d; g_tilde_2d = eta*H*gb_c_2d; phi_2d = phi*ones(1,Nx);

Boundary_data = fb - g_tilde + phi;
Full_data = f_2d - g_tilde_2d + phi_2d;

[val, error] = Hfun(Xp, 1, 0.0001);
real_soln = val/sqrt(3) - Xp;
real_soln_neg = real_soln(end:-1:1);
real_soln = [Xp;real_soln_neg];


figure(3)
set(gca,'fontsize',20);
plot3(zeros(N,1),Xp,fb_ini_s_data);
xlabel('x');ylabel('\mu');
title('given data f(x=0,\mu>0)'); xlim([0,0.5]); ylim([-1,1]);
print(gcf,'-depsc2','data/boundary_given_v.eps');
% pause;

figure(1)
set(gca,'fontsize',20);
plot(X,Full_data(:,1),'.',X,Boundary_data);
xlabel('\mu');ylabel('f_b');
title('f_b (v>0) = v, N = 80');
print(gcf,'-depsc2','data/boundary_v.eps');
% pause;

figure(2)
set(gca,'fontsize',20);
mesh(x_space,X,Full_data);
xlabel('x');ylabel('\mu');
title('f(x,v), N = 80'); xlim([0,0.5]);
print(gcf,'-depsc2','data/full_boundary_v.eps');
% pause;

save('data/exmaple_v');

% 
% figure(1)
% plot(X, fb, '.-', X, gb, 'ro-', X, phi, 'gx-');
% figure(2)
% set(gca,'fontsize',20);
% plot(X,Boundary_data,'o',X,real_soln);%Xp,fb_ini_s_data);
% % plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
% legend('numerical solution','exact solution','Location','NorthWest');
% xlabel('v');
% ylabel('f_b');
% title('f_b (v>0) = v, N = 80');
% print(gcf,'-depsc2','v.eps');