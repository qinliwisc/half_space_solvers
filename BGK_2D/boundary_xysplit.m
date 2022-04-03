%function boundary

%% preparation
Nx = 16; Ny = 16; u = 0; Nx_po = 16; Ny_po = 8;
indicator = 2;

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

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(Nx-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp - u; Xm = Xm - u; 

X = [Xp;Xm]; Wx = [Wp;Wm];
Hx = [Hp(:,1:Nx_po-1),Hp(:,1:Nx_po);...
    Hm(:,1:Nx_po-1),-Hm(:,1:Nx_po)];

vu = X + u;

[Y,Wy] = hermquad_symb(Ny-1);
[Y,ind] = sort(Y,'descend'); Wy = Wy(ind);
[Hy] = hermite_poly(Y, Ny-1) * (pi)^(1/4);
Wy = Wy/(Hy(:,1)'*diag(Wy)*Hy(:,1));

H_2d = kron(Hy(:,1:Ny_po),Hx);
W_2d = kron(Wy,Wx);
x_2d = kron(ones(length(Y),1),X);

[yy,xx] = meshgrid(Y,vu);

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
[fb_c] = BoundaryData_2d_xysplit(indicator, Nx, Ny, Nx_po, Ny_po, u, 0.1);
fb = H_2d*fb_c;
D = vchi'*diag(W_2d)*fb;

[C,g_c] = ConstructC_xysplit(Nx,Ny,Nx_po, Ny_po, u,0.1);

size_eta = length(index_p) + length(index_0);
eta = zeros(size_eta,1);

eta = C\D;


gb = H_2d*(g_c*eta);
phi = chi_data*eta;
boundary = fb - gb + phi;
boundary = reshape(boundary,length(X),length(Y));

exponential = exp(-xx.^2/2 - yy.^2/2);
Boundary = boundary.*exponential;

chi_s = chi(indicator);
fb_ini_half = feval(chi_s.fun,xx(1:length(Xp),:),yy(1:length(Xp),:));
fb_ini = feval(chi_s.fun,xx,yy);
fb_ini = reshape(fb_ini,length(X),length(Y));
Fb_ini = fb_ini.*exponential;

set(gca,'fontsize',20);
mesh(xx,yy,(Fb_ini - Boundary));
% legend('difference','location','southeast');
xlabel('x');
ylabel('y');
zlabel('f');zlim([-1 1]);
%title('boundary difference');
print(gcf,'-depsc2', ['data/indicator',num2str(indicator),'_diff.eps']);


set(gca,'fontsize',20);
mesh(xx,yy,Fb_ini);
% legend('given boundary','location','southeast');
xlabel('x');
ylabel('y');
zlabel('f');
%title('given boundary data');
print(gcf,'-depsc2', ['data/indicator',num2str(indicator),'_boundary.eps']);

set(gca,'fontsize',20);
mesh(xx,yy,Boundary);
% legend('recover','location','southeast');
xlabel('x');
ylabel('y');
zlabel('f');
%title('recovery');
print(gcf,'-depsc2', ['data/indicator',num2str(indicator),'_recover.eps']);

%return
