function boundary_xysplit_Maxwell(indicator,Nx_po,Ny_po,alpha1,alpha2,alpha3)

if (nargin == 0)
    alpha1 = 1; alpha2 = 0; alpha3 = 0;
    indicator = 4;
    Nx_po = 16; Ny_po = 8;
end
%% preparation
Nx = 16; Ny = 16; u = 0;

%% initial data
index_p = []; index_0 = [];

if (abs(u)<eps)
    index_p = [3];
    index_0 = [1,2];
elseif (u>0) 
    if (abs(u-sqrt(3/2))<eps)
        index_0 = [4];
        index_p = [1,2,3];
    elseif (u<sqrt(3/2))
        index_p = [1,2,3];
    elseif (u>sqrt(3/2))
        index_p = [1,2,3,4];
    end
elseif (u<0)
    if (abs(u+sqrt(3/2))<eps)
        index_0 = [3];
    elseif (u>-sqrt(3/2))
        index_p = [3];
    end
end

chi(1).fun = inline('(vx.^2+vy.^2-2)/sqrt(2*sqrt(pi))','vx','vy');
chi(2).fun = inline('sqrt(2)*vy/sqrt(sqrt(pi))','vx','vy');
chi(3).fun = inline('(vx.^2+vy.^2+2*vx)/2/sqrt(sqrt(pi))','vx','vy');
chi(4).fun = inline('(vx.^2+vy.^2-2*vx)/2/sqrt(sqrt(pi))','vx','vy');

fb_ini_given.fun = inline('(vx.^2+vy.^2-2)/sqrt(2*sqrt(pi))','vx','vy');

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(Nx-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

X = [Xp;Xm]; Wx = [Wp;Wm];
Hx = [Hp(:,1:Nx_po-1),Hp(:,1:Nx_po);...
    Hm(:,1:Nx_po-1),-Hm(:,1:Nx_po)];

[Y,Wy] = hermquad_symb(Ny-1);
[Y,ind] = sort(Y,'descend'); Wy = Wy(ind);
[Hy] = hermite_poly(Y, Ny-1) * (pi)^(1/4);
Wy = Wy/(Hy(:,1)'*diag(Wy)*Hy(:,1));

H_2d = kron(Hy(:,1:Ny_po),Hx);
W_2d = kron(Wy,Wx);
x_2d = kron(ones(length(Y),1),X);
xm_2d = kron(ones(size(Y)),Xm);
Wm_2d = kron(Wy,Wm);

[yy,xx] = meshgrid(Y,X);
[yy_pos,xx_pos] = meshgrid(Y,Xp);
[yy_neg_updown,xx_neg_updown] = meshgrid(Y,-Xp);
[yy_neg,xx_neg] = meshgrid(Y,Xm);

chip = []; chi0 = []; vchip = []; vchi0 = [];
for k=1:length(index_0)
    chi0_s = chi(index_0(k));
    chi0_d = feval(chi0_s.fun,xx,yy);
    chi0_d = reshape(chi0_d,length(x_2d),1);
    vchi0_d = x_2d.*chi0_d;

    vchi0 = horzcat(vchi0,vchi0_d);
    chi0 = horzcat(chi0,chi0_d);
end
for k=1:length(index_p)
    chip_s = chi(index_p(k));
    chip_d = feval(chip_s.fun,xx,yy);
    chip_d = reshape(chip_d,length(x_2d),1);
    vchip_d = x_2d.*chip_d;

    vchip = horzcat(vchip,vchip_d);
    chip = horzcat(chip,chip_d);
end

vchi = horzcat(vchip,vchi0);
chi_data = horzcat(chip,chi0);
[fb_c] = BoundaryData_2d_xysplit_Maxwell(fb_ini_given,indicator, Nx, ...
    Ny, Nx_po, Ny_po, u,alpha1,alpha2,alpha3, 0.1);
fb = H_2d*fb_c;
D = vchi'*diag(W_2d)*fb;

% [C,g_c] = ConstructC_xysplit(Nx,Ny,Nx_po, Ny_po, u,0.1);
[C,g_c] = ConstructC_xysplit_Maxwell(fb_ini_given,Nx,Ny,Nx_po, Ny_po,u, alpha1,alpha2,alpha3,0.1);

size_eta = length(index_p) + length(index_0);
eta = zeros(size_eta,1);

eta = C\D;

gb = H_2d*(g_c*eta);
phi = chi_data*eta;

boundary = fb - gb + phi;
boundary = reshape(boundary,length(X),length(Y));

exponential = exp(-xx.^2/2 - yy.^2/2);
Boundary = boundary.*exponential;

if indicator > 0
    chi_s = chi(indicator);
    fb_ini = feval(chi_s.fun,xx,yy);
    
	fb_ini1 = feval(chi_s.fun,xx_pos,yy_pos);
	fb_ini2 = feval(chi_s.fun,xx_neg_updown,yy_neg_updown);
	fb_ini3 = feval(chi_s.fun,xx_neg,yy_neg);
	fb_ini1 = reshape(fb_ini1',length(Xp)*length(Y),1);
	fb_ini2 = reshape(fb_ini2',length(Xm)*length(Y),1);
	fb_ini3 = reshape(fb_ini3',length(Xm)*length(Y),1);
	fb_ini3 = 2/sqrt(pi) * ones(length(Xp)*length(Y),1) * xm_2d'*diag(Wm_2d)*fb_ini3;
	fb_ini_half = fb_ini1 - alpha2 * fb_ini2 - alpha3 * fb_ini3;
    
    fb_ini = reshape(fb_ini,length(X),length(Y));
    Fb_ini = fb_ini.*exponential;
    fb_ini_half = reshape(fb_ini_half,length(Xp),length(Y));
    Fb_ini_half = fb_ini_half.*exponential(1:length(Xp),:);
elseif indicator == 0
    fb_ini = feval(fb_ini_given.fun,xx,yy);
    fb_ini_half = feval(fb_ini_given.fun,xx_pos,yy_pos);
    
    fb_ini = reshape(fb_ini,length(X),length(Y));
    Fb_ini = fb_ini.*exponential;
    fb_ini_half = reshape(fb_ini_half,length(Xp),length(Y));
    Fb_ini_half = fb_ini_half.*exponential(1:length(Xp),:);
end



set(gca,'fontsize',20);
mesh(xx_pos,yy_pos,Fb_ini_half);
% legend('Boundary difference','location','southeast');
xlabel('x','fontsize',20);
ylabel('y','fontsize',20);
zlabel('h','fontsize',20);
saveas(gcf,['data/maxwell/indicator',num2str(indicator),'_Dirichlet.eps'],'epsc');
%title('boundary given Dirichlet');
% print(gcf,'-depsc2', ['data/maxwell/indicator',num2str(indicator),'_Dirichlet.eps']);
% print(gcf,'-dpdf', ['data/maxwell/indicator',num2str(indicator),'_Dirichlet.pdf']);

set(gca,'fontsize',20);
mesh(xx,yy,(-Fb_ini + Boundary));
% legend('Boundary difference','location','southeast');
xlabel('x','fontsize',20);
ylabel('y','fontsize',20);
zlabel('f_h-X_-','fontsize',20); zlim([-1 1]);
saveas(gcf,['data/maxwell/indicator',num2str(indicator),'_diff.eps'],'epsc');
%title('boundary difference');
% print(gcf,'-depsc2', ['data/maxwell/indicator',num2str(indicator),'_diff.eps']);
% print(gcf,'-dpdf', ['data/maxwell/indicator',num2str(indicator),'_diff.pdf']);

if indicator > 0
set(gca,'fontsize',20);
mesh(xx,yy,Fb_ini);
% legend('Boundary given','location','southeast');
xlabel('x','fontsize',20);
ylabel('y','fontsize',20);
zlabel('X_-','fontsize',20);
%title('analytical solution');
% print(gcf,'-deps2', ['data/maxwell/indicator',num2str(indicator),'_boundary.eps']);
% print(gcf,'-dpdf', ['data/maxwell/indicator',num2str(indicator),'_boundary.pdf']);
saveas(gcf,['data/maxwell/indicator',num2str(indicator),'_boundary.eps'],'epsc');
end

set(gca,'fontsize',20);
mesh(xx,yy,Boundary);
% legend('recover','location','southeast');
xlabel('x','fontsize',20);
ylabel('y','fontsize',20);
zlabel('f_h','fontsize',20);
%title('recovery');
% print(gcf,'-depsc2', ['data/maxwell/indicator',num2str(indicator),'_recover.eps']);
% print(gcf,'-dpdf', ['data/maxwell/indicator',num2str(indicator),'_recover.pdf']);
saveas(gcf,['data/maxwell/indicator',num2str(indicator),'_recover.eps'],'epsc');

return
