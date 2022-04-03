function [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = ...
    boundary_Maxwell(indicator,N,alpha1,alpha2,alpha3,fb_ini)

%% preparation
if (nargin == 0)
    clear;
    N = 16; N_omega = 8; indicator = 1;
    alpha1 = 0.3; alpha2 = 0.3; alpha3 = 0.4;
    fb_ini.fun = inline('v - v + omega.^(-1/2)','v','omega');
end


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
alpha = omega.^(1/2); alpha = alpha(:);
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

%% discretization on v
[omegaomega,xx] = meshgrid(omega,X);

%% recovry
% [fb_c] = BoundaryData_Maxwell(N,N_omega,indicator,alpha1,alpha2,alpha3,fb_ini);
% [gb_c] = BoundaryData_Maxwell(N,N_omega,1,alpha1,alpha2,alpha3,fb_ini);
[fb_c,h] = BoundaryData_Maxwell_pointwise(N,N_omega,indicator,alpha1,alpha2,alpha3,fb_ini);
[gb_c,~] = BoundaryData_Maxwell_pointwise(N,N_omega,1,alpha1,alpha2,alpha3,fb_ini);

fb = H*fb_c; gb = H*gb_c;

C = X_2d'*diag(Weight)*gb;
D = X_2d'*diag(Weight)*fb;
eta = D/C;

g_tilde = gb*eta; phi = chi_data*eta;

Boundary_data = fb - g_tilde + phi;

fb = reshape(fb,length(X),N_omega);
g_tilde = reshape(g_tilde,length(X),N_omega);
phi = reshape(phi,length(X),N_omega);
Boundary_data = reshape(Boundary_data,length(X),N_omega);
h = reshape(h,length(Xp),N_omega);

if indicator > 0
    fb_ini_data = ones(length(X),1)*alpha_inv';
elseif indicator == 0
    fb_ini_data = feval(fb_ini.fun,xx,omegaomega);
end

figure(2)
set(gca,'fontsize',20);
mesh(omega,Xp,h);
xlabel('\omega','fontsize',20);
ylabel('\mu','fontsize',20); %zlabel('h');
%title('recovery');
saveas(gcf,['data/maxwell_pointwise/indicator',num2str(indicator),'_h.eps'],'epsc');
% print(gcf,'-depsc2',['data/maxwell_pointwise/indicator',num2str(indicator),'_h.eps']);pause;
% print(gcf,'-dpdf',['data/maxwell_pointwise/indicator',num2str(indicator),'_h.pdf']);pause;
% close(figure(2))

figure(3)
set(gca,'fontsize',20);
% gca.FontSize = 20;
mesh(omega,X,Boundary_data);
xlabel('\omega','fontsize',20);
ylabel('\mu','fontsize',20); %zlabel('f_h');
saveas(gcf,['data/maxwell_pointwise/indicator',num2str(indicator),'_recover.eps'],'epsc')
%title('recovery');
% print(gcf,'-depsc2',['data/maxwell_pointwise/indicator',num2str(indicator),'_recover.eps']);pause;
% print(gcf,'-dpdf',['data/maxwell_pointwise/indicator',num2str(indicator),'_recover.pdf']);pause;
% close(figure(3))

figure(4)
set(gca,'fontsize',20);
mesh(omega,X,(-fb_ini_data + Boundary_data));
xlabel('\omega','fontsize',20);
ylabel('\mu','fontsize',20);
saveas(gcf,['data/maxwell_pointwise/indicator',num2str(indicator),'_diff.eps'],'epsc')
% zlabel('f_h-h');
%zlim([-1 1]);
%title('boundary difference');
% print(gcf,'-depsc2',['data/maxwell_pointwise/indicator',num2str(indicator),'_recover.eps']);pause;
% print(gcf,'-dpdf', ['data/maxwell_pointwise/indicator',num2str(indicator),'_diff.pdf']);pause;
% close(figure(4))

figure(5)
set(gca,'fontsize',20);
mesh(omega,X,fb_ini_data);
xlabel('\omega','fontsize',20);
ylabel('\mu','fontsize',20);
% zlabel('h','fontsize',20);
saveas(gcf, ['data/maxwell_pointwise/indicator',num2str(indicator),'_boundary.eps'],'epsc');
%title('given boundary data');
% print(gcf,'-dpdf',['data/maxwell_pointwise/indicator',num2str(indicator),'_boundary.pdf']);pause;
% print(gcf,'-depsc2', ['data/maxwell_pointwise/indicator',num2str(indicator),'_boundary.eps']);pause;
% close(figure(5))