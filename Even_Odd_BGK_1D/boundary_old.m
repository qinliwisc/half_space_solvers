function [x_dense,Fb_ini_dense,Phi_dense,Boundary_data_dense]=boundary(N,u,fb_ini_s)
% function [Fb_ini,Phi,Boundary_data,vu,X,Xp]=boundary(N,u,fb_ini_s)

%% preparation
if (nargin == 0)
    clear;
    N = 16; u = -sqrt(1.5);%0.5;
    fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
end

%% initial data
index_p = []; index_0 = [];

if (abs(u)<eps)
    index_p = [2];
    index_0 = [1];
elseif (u>0) 
    if (abs(u-sqrt(3/2))<eps)
        index_0 = [3];
        index_p = [1,2];
    elseif (u<sqrt(3/2))
        index_p = [1,2];
    elseif (u>sqrt(3/2))
        index_p = [1,2,3];
    end
elseif (u<0)
    if (abs(u+sqrt(3/2))<eps)
        index_0 = [2];
    elseif (u>-sqrt(3/2))
        index_p = [2];
    elseif (u<-sqrt(3/2))
        error('NO solution in this case');
        return
    end
end

chi(1).fun = inline('(2*v.^2-3)/sqrt(6)');
chi(2).fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
chi(3).fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
% fb_ini_s.fun = chi(1).fun;

%% discretization on v: inner product quadrature
[Xp, Wp] = half_hermquad_shift_symb(N-1, u/2);
Wp = exp(-u^2/4) / sqrt(pi) * Wp;
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 

dense_n = 100; dx = (Xp(1)-Xp(end))/dense_n;
x_dense_p = Xp(end) + [0:dense_n]*dx; x_dense_p = x_dense_p(:);
x_dense_p = flipud(x_dense_p);
[Hp_dense] = half_hermite_shift_poly(x_dense_p, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 

Xp = Xp - u; % shift to center at -u;
x_dense_p = x_dense_p-u;

[Xm, Wm] = half_hermquad_shift_symb(N-1, -u/2);
Wm = exp(-u^2/4) / sqrt(pi) * Wm;
[Xm,ind] = sort(Xm,'ascend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(Xm, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2);

dense_n = 100; dx = (Xm(end)-Xm(1))/dense_n;
x_dense_m = Xm(1) + [0:dense_n]*dx; x_dense_m = x_dense_m(:);
% x_dense_m = flipud(x_dense_m);
[Hm_dense] = half_hermite_shift_poly(x_dense_m, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 

Xm = - Xm - u; % shift to center at -u;
x_dense_m = -x_dense_m - u;

X = [Xp; Xm]; x_dense = [x_dense_p;x_dense_m];
W = [Wp; Wm]; 


H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
H_dense = [Hp_dense(:,1:end-1),Hp_dense;Hm_dense(:,1:end-1),-Hm_dense];

vu = X + u;

%% recovery
chip = []; chi0 = [];
chip_dense = []; chi0_dense = [];
vuchip = []; vuchi0 = [];
for k=1:length(index_0)
    chi0_s = chi(index_0(k));
    vuchi0_d = vu.*feval(chi0_s.fun,X); vuchi0_d = vuchi0_d(:);
    vuchi0 = horzcat(vuchi0,vuchi0_d);
    chi0_d = feval(chi0_s.fun,X); chi0_d = chi0_d(:);
    chi0 = horzcat(chi0,chi0_d);
    chi0_d_dense = feval(chi0_s.fun,x_dense); chi0_d_dense = chi0_d_dense(:);
    chi0_dense = horzcat(chi0_dense,chi0_d_dense);
end
for k=1:length(index_p)
    chip_s = chi(index_p(k));
    vuchip_d = vu.*feval(chip_s.fun,X); vuchip_d = vuchip_d(:);
    vuchip = horzcat(vuchip,vuchip_d);
    chip_d = feval(chip_s.fun,X); chip_d = chip_d(:);
    chip = horzcat(chip,chip_d);
    chip_d_dense = feval(chip_s.fun,x_dense); chip_d_dense = chip_d_dense(:);
    chip_dense = horzcat(chip_dense,chip_d_dense);
end
vuchi = horzcat(vuchip,vuchi0);
chi_data = horzcat(chip,chi0);
chi_data_dense = horzcat(chip_dense,chi0_dense);

[fb_c]=BoundaryData_c(N,u,0,X,W,H,fb_ini_s); fb = H*fb_c;
Fb = fb/(pi^(1/4)).*exp(-vu.^2/2);
D = vuchi'*diag(W)*fb;

[C,g_c_d] = ConstructC(N,u,index_p,index_0,X,W,H,fb_ini_s);
eta = C\D;

fb_ini_d = feval(fb_ini_s.fun,Xp);
Fb_ini = fb_ini_d/(pi^(1/4)).*exp(-(Xp).^2/2);
phi = chi_data*eta;
Phi = phi/(pi^(1/4)).*exp(-X.^2/2); 

% plot(Xp,Fb_ini,'.-.',vu,Phi,'o-');
% legend('real: u mode','f_\infty','location','southeast');
% xlabel('v');
% ylabel('f_b/f_\infty');
% title('given f_b (v>0) = u mode');

g_tilde = H*(g_c_d*eta);
G_tilde = g_tilde/(pi^(1/4)).*exp(-vu.^2/2);
Boundary_data = Fb - G_tilde + Phi;
% Boundary_data = boundary_data/(pi^(1/4)).*exp(-vu.^2/2);

% plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
% legend('real: u mode','recover, full','location','southeast');
% xlabel('v');
% ylabel('f_b');
% title('given f_b (v>0) = u mode');

fb_ini_d_dense = feval(fb_ini_s.fun,x_dense);
Fb_ini_dense = fb_ini_d_dense/(pi^(1/4)).*exp(-(x_dense).^2/2);
phi_dense = chi_data_dense*eta;
Phi_dense = phi_dense/(pi^(1/4)).*exp(-(x_dense).^2/2);
fb_dense = H_dense*fb_c;
Fb_dense = fb_dense/(pi^(1/4)).*exp(-(x_dense+u).^2/2);
g_tilde_dense = H_dense*(g_c_d*eta);
G_tilde_dense = g_tilde_dense/(pi^(1/4)).*exp(-(x_dense+u).^2/2);
Boundary_data_dense = Fb_dense - G_tilde_dense + Phi_dense;

% plot(x_dense,Fb_ini_dense,'.-.',x_dense,Phi_dense,'r',vu,Boundary_data,'o-');
% legend('real: u mode','f at infinity','recover, full','location','southeast');
% xlabel('v');
% ylabel('f(v)');
% title('given f_b (v>0) = u mode');

return