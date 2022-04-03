% function [fb_c,g_coeff,eta,index_total] = boundary_coeff(N,u,fb_ini_s)
function [fb_c,g_coeff,eta,phi_coeff,coeff,X,Boundary_data] = boundary_coeff(N,u,fb_ini_s)

%% preparation
if (nargin == 0)
    clear;
    N = 36; 
    u = 0;
    fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');
	fb_ini_s.fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');
	fb_ini_s.fun = inline('v.^3');
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
index_total = [index_p,index_0];

chi(1).fun = inline('(2*v.^2-3)/sqrt(6)');
chi(2).fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
chi(3).fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');

filenamex = 'prepare/x_';filenamex = strcat(filenamex,num2str(N));
filenamex = strcat(filenamex,'.dat');
filenamew = 'prepare/w_';filenamew = strcat(filenamew,num2str(N));
filenamew = strcat(filenamew,'.dat');

%% discretization on v: inner product quadrature
if u==0
    [Xp] = load(filenamex);Wp = load(filenamew);
else
    [Xp, Wp] = half_hermquad_shift_symb(N-1, u/2);
end
Wp = exp(-u^2/4) / sqrt(pi) * Wp;
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2); 
Xp = Xp - u; % shift to center at -u;

if u==0
    [Xm] = load(filenamex);Wm = load(filenamew);
else
    [Xm, Wm] = half_hermquad_shift_symb(N-1, -u/2);
end
Wm = exp(-u^2/4) / sqrt(pi) * Wm;
[Xm,ind] = sort(Xm,'ascend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(Xm, [0:1:N-1], 0)*(pi)^(1/4)/sqrt(2);
Xm = - Xm - u; % shift to center at -u;

X = [Xp; Xm]; W = [Wp; Wm]; 

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

vu = X + u;

%% recovery
chi_data = []; vuchi = [];
for k=1:length(index_total)
    chi_s = chi(index_total(k));
    vuchi_d = vu.*feval(chi_s.fun,X); vuchi_d = vuchi_d(:);
    vuchi = horzcat(vuchi,vuchi_d);
    chi_d = feval(chi_s.fun,X); chi_d = chi_d(:);
    chi_data = horzcat(chi_data,chi_d);
end

[fb_c] = BoundaryData_c(N,u,0,X,W,H,fb_ini_s); fb = H*fb_c;
Fb = fb/(pi^(1/4)).*exp(-vu.^2/2);
D = vuchi'*diag(W)*fb;

[C,g_c_d] = ConstructC(N,u,index_p,index_0,X,W,H,fb_ini_s);
eta = C\D; g_coeff = g_c_d*eta;

% fb_ini_d = feval(fb_ini_s.fun,Xp);
% Fb_ini = fb_ini_d/(pi^(1/4)).*exp(-(Xp).^2/2);

phi = chi_data*eta;
phi_coeff = H'*diag(W)*phi;
Phi = phi/(pi^(1/4)).*exp(-X.^2/2); 

g_tilde = H*(g_coeff);
G_tilde = g_tilde/(pi^(1/4)).*exp(-vu.^2/2);
Boundary_data = Fb - G_tilde + Phi;

coeff = fb_c - g_coeff + phi_coeff;
plot(X,G_tilde,'.-.',X,Fb,'o');
% plot(Xp,Fb_ini,'.-.',X,Boundary_data,'o-');
legend('g','f');
xlabel('v');
ylabel('f_b');
title('given f_b (v>0) = u mode');

% fb_ini_d_dense = feval(fb_ini_s.fun,x_dense);
% Fb_ini_dense = fb_ini_d_dense/(pi^(1/4)).*exp(-(x_dense).^2/2);
% phi_dense = chi_data_dense*eta;
% Phi_dense = phi_dense/(pi^(1/4)).*exp(-(x_dense).^2/2);
% fb_dense = H_dense*fb_c;
% Fb_dense = fb_dense/(pi^(1/4)).*exp(-(x_dense+u).^2/2);
% g_tilde_dense = H_dense*(g_c_d*eta);
% G_tilde_dense = g_tilde_dense/(pi^(1/4)).*exp(-(x_dense+u).^2/2);
% Boundary_data_dense = Fb_dense - G_tilde_dense + Phi_dense;

return