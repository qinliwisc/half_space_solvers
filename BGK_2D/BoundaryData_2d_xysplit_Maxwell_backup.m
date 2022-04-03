function [c]=BoundaryData_2d_xysplit_Maxwell(index, Nx, Ny, Nx_po, Ny_po, u, alpha1,alpha2,alpha3,dampcoef)

if (nargin == 0)
    Nx = 32; Ny = 16; u = 0; dampcoef = 0.1; index = 1;
    Nx_po = 16; Ny_po = 8;
    alpha1 = 0.3; alpha2 = 0.3; alpha3 = 0.4;
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
Hm_updown = flipud(Hm);

X = [Xp;Xm]; Wx = [Wp;Wm];
Hx = [Hp(:,1:Nx_po - 1),Hp(:,1:Nx_po);...
    Hm(:,1:Nx_po - 1),-Hm(:,1:Nx_po)];


[Y,Wy] = hermquad_symb(Ny-1);
[Y,ind] = sort(Y,'descend'); Wy = Wy(ind);
[Hy] = hermite_poly(Y, Ny-1) * (pi)^(1/4);
Wy = Wy/(Hy(:,1)'*diag(Wy)*Hy(:,1));

[yy,xx] = meshgrid(Y,X);


%% polynomials basis
H_2d = kron(Hy(:,1:Ny_po),Hx);
x_2d = kron(ones(size(Y)),Xp);
Wp_2d = kron(Wy,Wp);
xm_2d = kron(ones(size(Y)),Xm);
Wm_2d = kron(Wy,Wm);
Weight = kron(Wy,Wx);

Hpp = kron(Hy(:,1:Ny_po),Hp(:,1:Nx_po-1));
H_pos = kron(Hy(:,1:Ny_po),[Hp(:,1:Nx_po-1),Hp(:,1:Nx_po)]);
H_neg_updown = kron(Hy(:,1:Ny_po),[Hm_updown(:,1:Nx_po-1),Hm_updown(:,1:Nx_po)]);
H_neg = kron(Hy(:,1:Ny_po),[Hm(:,1:Nx_po-1),Hm(:,1:Nx_po)]);

%% initial data
if index == 1
    fb_ini = chi_0x(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);

    fb_ini = chi_0x(xx,yy);
    fb_ini = reshape(fb_ini,length(X)*length(Y),1);
    c_ini = H_2d'*diag(Weight)*fb_ini;
    fb_ini = (H_pos - alpha2 * H_neg_updown)* c_ini;
    fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
        ones(length(x_2d),1) * xm_2d'*diag(Wm_2d)*H_neg * c_ini;
    fb_ini = fb_ini/alpha1;
    
%     fb_ini_updown = chi_0x(-xx(1:length(Xp),:),yy(1:length(Xp),:));
%     fb_ini_updown = reshape(fb_ini_updown,length(Xp)*length(Y),1);
%     
%     fb_ini_neg = chi_0x(xx(1+length(Xp):end,:),yy(1+length(Xp):end,:));
%     fb_ini_neg = reshape(fb_ini_neg,length(Xm)*length(Y),1);
%     
%     fb_ini = fb_ini - alpha2 * fb_ini_updown;
%     fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
%         ones(length(x_2d),1) * (xm_2d'*diag(Wm_2d)*fb_ini_neg);
%     fb_ini = fb_ini/alpha1;
elseif index == 2
    fb_ini = chi_0y(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
 
    fb_ini = chi_0y(xx,yy);
    fb_ini = reshape(fb_ini,length(X)*length(Y),1);
    c_ini = H_2d'*diag(Weight)*fb_ini;
    fb_ini = (H_pos - alpha2 * H_neg_updown)* c_ini;
    fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
        ones(length(x_2d),1) * xm_2d'*diag(Wm_2d)*H_neg * c_ini;
    fb_ini = fb_ini/alpha1;

%     fb_ini_updown = chi_0y(-xx(1:length(Xp),:),yy(1:length(Xp),:));
%     fb_ini_updown = reshape(fb_ini_updown,length(Xp)*length(Y),1);
%     
%     fb_ini_neg = chi_0y(xx(1+length(Xp):end,:),yy(1+length(Xp):end,:));
%     fb_ini_neg = reshape(fb_ini_neg,length(Xm)*length(Y),1);
%     
%     fb_ini = fb_ini - alpha2 * fb_ini_updown;
%     fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
%         ones(length(x_2d),1) * (xm_2d'*diag(Wm_2d)*fb_ini_neg);
%     fb_ini = fb_ini/alpha1;
elseif index == 3
    fb_ini = chi_p(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    
%     fb_ini_updown = chi_p(-xx(1:length(Xp),:),yy(1:length(Xp),:));
%     fb_ini_updown = reshape(fb_ini_updown,length(Xp)*length(Y),1);
%     
%     fb_ini_neg = chi_p(xx(1+length(Xp):end,:),yy(1+length(Xp):end,:));
%     fb_ini_neg = reshape(fb_ini_neg,length(Xm)*length(Y),1);
%     
%     fb_ini = fb_ini - alpha2 * fb_ini_updown;
%     fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
%         ones(length(x_2d),1) * (xm_2d'*diag(Wm_2d)*fb_ini_neg);
%     fb_ini = fb_ini/alpha1;
    
    fb_ini = chi_p(xx,yy);
    fb_ini = reshape(fb_ini,length(X)*length(Y),1);
    c_ini = H_2d'*diag(Weight)*fb_ini;
    fb_ini = (H_pos - alpha2 * H_neg_updown)* c_ini;
    fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
        ones(length(x_2d),1) * xm_2d'*diag(Wm_2d)*H_neg * c_ini;
    fb_ini = fb_ini/alpha1;

elseif index == 4
    fb_ini = chi_m(xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini = reshape(fb_ini,length(Xp)*length(Y),1);
    
    fb_ini_updown = chi_m(-xx(1:length(Xp),:),yy(1:length(Xp),:));
    fb_ini_updown = reshape(fb_ini_updown,length(Xp)*length(Y),1);
    
    fb_ini_neg = chi_m(xx(1+length(Xp):end,:),yy(1+length(Xp):end,:));
    fb_ini_neg = reshape(fb_ini_neg,length(Xm)*length(Y),1);
    
    fb_ini = fb_ini - alpha2 * fb_ini_updown;
    fb_ini = fb_ini + 2/sqrt(pi) * alpha3 *...
        ones(length(x_2d),1) * (xm_2d'*diag(Wm_2d)*fb_ini_neg);
    fb_ini = fb_ini/alpha1;
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
term1 = Hpp' * diag(Wp_2d.*x_2d) * H_pos;
term1 = term1 - alpha2 * Hpp' * diag(Wp_2d.*x_2d) * H_neg_updown;
term1 = term1 + 2/sqrt(pi) * alpha3 * Hpp' * diag(Wp_2d.*x_2d) * ones(length(x_2d),1) * xm_2d'*diag(Wm_2d)*H_neg;
A(Ny_po*Nx_po+1:end,:) = term1;
b(Ny_po*Nx_po+1:end, :) = alpha1 * Hpp' * diag(Wp_2d.*x_2d) * fb_ini;

c = A \ b; 

return
