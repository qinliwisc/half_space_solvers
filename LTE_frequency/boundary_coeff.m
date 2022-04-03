% function [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = ...
%     boundary_coeff(N,fb_ini)
% 
% %% preparation
% if (nargin == 0)
    clear;
    N = 16; N_omega = 8; index = 1;
    if index == 1
        fb_ini.fun = inline('v - v + omega.^(-1/2)','v','omega');
    elseif index == 2
        fb_ini.fun = inline('v - v + omega.^(1/2)','v','omega');
    end
% end


%% discretization on v: inner product quadrature
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
Hp = half_legendre_poly(Xp,[0:1:N-1]); Hp = Hp(ind,:);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

X = [Xp; Xm]; W = [Wp; Wm]; W = W/2;

H_X = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

omega = [1:1:N_omega];

[omegaomega,xx] = meshgrid(omega,X);

alpha = omega.^(1/2);alpha = alpha(:);
beta = omega.*exp(-omega/1000); beta = beta(:);

beta_inv = ones(size(beta))./beta;
alpha_inv = ones(size(alpha))./alpha;
H = kron(diag(sqrt(beta_inv)),H_X);
Weight = kron(beta,W); Weight = Weight / sum(Weight);
iii = H'*diag(Weight)*H;
H = H/sqrt(iii(1,1));
X_2d = kron(ones(N_omega,1),X);

chi_data = ones(length(X),1)*alpha_inv';
chi_data = reshape(chi_data,N_omega*length(X),1);

%% recovry
[fb_c] = BoundaryData(N,N_omega,0,fb_ini);
[gb_c] = BoundaryData(N,N_omega,1,fb_ini);

fb = H*fb_c; gb = H*gb_c;

C = X_2d'*diag(Weight)*gb;
D = X_2d'*diag(Weight)*fb;
eta = D/C;

g_tilde = gb*eta;
phi = chi_data*eta;

Boundary_data = fb - g_tilde + phi;

fb = reshape(fb,length(X),N_omega);
g_tilde = reshape(g_tilde,length(X),N_omega);
phi = reshape(phi,length(X),N_omega);
Boundary_data = reshape(Boundary_data,length(X),N_omega);

fb_ini_given = feval(fb_ini.fun,xx,omegaomega);

% figure(1)
% set(gca,'fontsize',20);
% mesh(omega,X,fb);
% title(['N = ' num2str(2*N-1) ', N_\omega = 8']);
% print(gcf,'-depsc2',['data_galerkin/f',num2str(index),'_damped.eps']);
% close(figure(1))
% 
% figure(2)
% set(gca,'fontsize',20);
% mesh(omega,X,g_tilde);
% title(['N = ' num2str(2*N-1) ', N_\omega = 8']);
% print(gcf,'-depsc2',['data_galerkin/g',num2str(index),'_damped.eps']);
% close(figure(2))
% 
figure(3)
set(gca,'fontsize',20);
mesh(omega,X,Boundary_data);
% title(['recovery']);
xlabel('\omega','fontsize',20);ylabel('\mu','fontsize',20);
% zlabel('f_h','fontsize',20);
saveas(gcf,['data/recover',num2str(index),'.eps'],'epsc');
% print(gcf,'-dpdf',['data/recover',num2str(index),'.pdf']);
close(figure(3))

figure(4)
set(gca,'fontsize',20);
mesh(omega,X,Boundary_data - fb_ini_given);
% title(['boundary difference']);
xlabel('\omega','fontsize',20);
ylabel('\mu','fontsize',20);
%zlabel('f_h-h','fontsize',20);
zlim([-1 1]);
% caxis([-1 1]);
saveas(gcf,['data/diff',num2str(index),'.eps'],'epsc');
% print(gcf,'-dpdf',['data/diff',num2str(index),'.pdf']);
close(figure(4))