function [c]=BoundaryData_2d_xysplit(index, Nx, Ny, Nx_po, Ny_po, u, dampcoef)

if (nargin == 0)
    Nx = 32; Ny = 16; u = 0; dampcoef = 0.1; index = 1;
    Nx_po = 16; Ny_po = 8;
end

% 2D weight = exp(-(vx^2+vy^2)/2)
chi_0x = inline('(vx.^2+vy.^2-2)/sqrt(2*sqrt(pi))','vx','vy');
chi_0y = inline('sqrt(2)*vy/sqrt(sqrt(pi))','vx','vy');
chi_p = inline('(vx.^2+vy.^2+2*vx)/2/sqrt(sqrt(pi))','vx','vy');
chi_m = inline('(vx.^2+vy.^2-2*vx)/2/sqrt(sqrt(pi))','vx','vy');

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(Nx-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp-u; Xm = Xm-u; 

X = [Xp;Xm]; Wx = [Wp;Wm];
Hx = [Hp(:,1:Nx_po - 1),Hp(:,1:Nx_po);...
    Hm(:,1:Nx_po - 1),-Hm(:,1:Nx_po)];

vu = X + u;

[Y,Wy] = hermquad_symb(Ny-1);
[Y,ind] = sort(Y,'descend'); Wy = Wy(ind);
[Hy] = hermite_poly(Y, Ny-1) * (pi)^(1/4);
Wy = Wy/(Hy(:,1)'*diag(Wy)*Hy(:,1));

[yy,xx] = meshgrid(Y,vu);


%% polynomials basis
H_2d = kron(Hy(:,1:Ny_po),Hx);
Wp_2d = kron(Wy,Wp);
x_2d = kron(ones(size(Y)),vu(1:length(Xp)));

Hpp = kron(Hy(:,1:Ny_po),Hp(:,1:Nx_po-1));
H_pos = kron(Hy(:,1:Ny_po),[Hp(:,1:Nx_po-1),Hp(:,1:Nx_po)]);


%% initial data
if index == 1
    fb_ini = chi_0x(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    fb_ini_true = chi_0x(xx,yy);
elseif index == 2
    fb_ini = chi_0y(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    fb_ini_true = chi_0y(xx,yy);
elseif index == 3
    fb_ini = chi_p(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    fb_ini_true = chi_p(xx,yy);
elseif index == 4
    fb_ini = chi_m(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    fb_ini_true = chi_m(xx,yy);
end

%% linear system
% generalized eigenvalue
[V,D,Amat,Bmat] = GeneralizedEigen_2d_xysplit(...
    Nx,Ny,u,dampcoef,Nx_po,Ny_po);
if (sum(D>-1e-13) ~= Ny_po*Nx_po)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

% consistency check for the eignevectors
if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
    error('Something wrong in the generalized eigenvalue');
end

% full system
A = zeros((2*Nx_po-1)*Ny_po,(2*Nx_po-1)*Ny_po);
b = zeros((2*Nx_po-1)*Ny_po, 1);

% The first N(2N-1) rows are constraints from the growing and zero modes
A(1:Ny_po*Nx_po, :) = Vcons.'*Bmat; 
b(1:Ny_po*Nx_po, :) = 0;

%% Collocation
A(Ny_po*Nx_po+1:end,:) = Hpp' * diag(Wp_2d.*x_2d) * H_pos;
b(Ny_po*Nx_po+1:end, :) = Hpp' * diag(Wp_2d.*x_2d) * fb_ini;

c = A \ b; 

exponential = exp(-xx.^2/2 - yy.^2/2);


fb = H_2d * c;
fb = reshape(fb,length(X),length(Y));
Fb = fb.*exponential;

fb_ini = reshape(fb_ini,length(Xp),length(Y));
Fb_ini = fb_ini.*exponential(1:length(Xp),:);
Fb_ini_true = fb_ini_true.*exponential;

diff = Fb - Fb_ini_true;

% figure(1)
% mesh(xx,yy,Fb);
% xlabel('vx');
% ylabel('vy');
% zlabel('f');
% figure(2)
% mesh(xx(1:length(Xp),:),yy(1:length(Xp),:),diff./Fb_ini)
% legend('difference','location','southeast');
% xlabel('vx');
% ylabel('vy');
% zlabel('f');
% title('given f_b (vx>0) = 1 mode');

return
