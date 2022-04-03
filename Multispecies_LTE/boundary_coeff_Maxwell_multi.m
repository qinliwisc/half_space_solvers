function [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = ...
    boundary_coeff_Maxwell_multi(N,fb_ini_s1,fb_ini_s2,sigma,alpha,casenum)

%% preparation
if (nargin == 0)
    clear;
    N = 80;
    fb_ini_s1.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s1.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
	fb_ini_s1.fun = inline('v');
    
    fb_ini_s2.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s2.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
	fb_ini_s2.fun = inline('v');
    
    sigma11 = 0.5; sigma12 = 0.5; sigma21 = 0.5; sigma22 = 0.5;
%     sigma11 = 0.25; sigma12 = sqrt(0.25*0.75); sigma21 = sqrt(0.25*0.75); sigma22 = 0.75;
%     sigma11 = 0.25; sigma12 = 0.25; sigma21 = 0.75; sigma22 = 0.75;
%     sigma11 = 1; sigma12 = 0; sigma21 = 0; sigma22 = 1;
    sigma = [sigma11,sigma12;sigma21,sigma22];
    
%      alpha = [1,0,0];
%     alpha = [0.9,0,0.1];
    alpha = [0.1,0.45,0.45];
%     alpha = [0.8,0.1,0.1];
end


%% discretization on v: inner product quadrature
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_legendre_poly(Xp,[0:1:N-1]);

fb_ini_s_data = [feval(fb_ini_s1.fun,Xp),feval(fb_ini_s2.fun,Xp)];

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

X = [Xp; Xm]; W = [Wp; Wm]; W = W/2;

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
H_tilde = [H,zeros(size(H));zeros(size(H)),H];

%% recovry
chi_data = [sigma(1,2)*ones(length(X),1);(1-sigma(1,1))*ones(length(X),1)];

[fb_c] = BoundaryData_Maxwell_multi(N,0,fb_ini_s1,fb_ini_s2,alpha,sigma);
[gb_c] = BoundaryData_Maxwell_multi(N,1,fb_ini_s1,fb_ini_s2,alpha,sigma);

fb = H_tilde*fb_c;% plot(X,fb);pause;
gb = H_tilde*gb_c;% plot(X,gb);pause;

C = [sigma(1,2)*X;(1-sigma(1,1))*X]'*diag([W;W])*gb;
D = [sigma(1,2)*X;(1-sigma(1,1))*X]'*diag([W;W])*fb;

eta = D/C;

g_tilde = gb*eta;
phi = chi_data*eta*alpha(1)/(1 - alpha(2) - alpha(3));

Boundary_data = fb - g_tilde + phi;

fb = [fb(1:length(X)),fb(length(X)+1:end)];
gb = [gb(1:length(X)),gb(length(X)+1:end)];
g_tilde = [g_tilde(1:length(X)),g_tilde(length(X)+1:end)];
phi = [phi(1:length(X)),phi(length(X)+1:end)];
Boundary_data = [Boundary_data(1:length(X)),Boundary_data(length(X)+1:end)];


figure(1)
plot(X, fb(:,1), '.-', X, gb(:,1), 'ro-', X, phi(:,1), 'gx-');
pause;
close(figure(1))

figure(2)
plot(X, fb(:,2), '.-', X, gb(:,2), 'ro-', X, phi(:,2), 'gx-');
pause;
close(figure(2))

figure(3)
set(gca,'fontsize',20);
plot(X,Boundary_data(:,1),'ro-', Xp,fb_ini_s_data(:,1),...
    X,Boundary_data(:,2),'gx-', Xp,fb_ini_s_data(:,2));
% plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
legend('solution 1','boundary data 1','solution 2','boundary data 2','Location','NorthWest');
xlabel('v');
ylabel('f_b');
title('f_{b1} (v>0) = v, f_{b2} (v>0) = v, N = 80');
print(gcf,'-depsc2',['test_v_2v_',num2str(casenum),'.eps']);
pause;
close(figure(3))

% figure(4)
% set(gca,'fontsize',20);
% plot(X,Boundary_data(:,2),'o',Xp,fb_ini_s_data(:,2));
% % plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
% legend('solution','boundary data','Location','NorthWest');
% xlabel('v');
% ylabel('f_b');
% title('f_b (v>0) = v, N = 80');
% print(gcf,'-depsc2','v_2.eps');