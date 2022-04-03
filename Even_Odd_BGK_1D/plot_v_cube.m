clear;

fb_ini_s_0.fun = inline('v.^3');

N = 36;
u = 0;

[x_dense_0,Fb_ini_dense_0,Phi_dense_0,Boundary_data_0]=boundary(N,u,fb_ini_s_0,8);

x_dense_half_0 = x_dense_0(x_dense_0>0);
Fb_ini_dense_half_0 = Fb_ini_dense_0(x_dense_0>0);

filename = 'data/Filtering2/u_convergence';
filename_0 = strcat(filename,'/zero_N36.mat');
save(filename_0,'x_dense_0','Fb_ini_dense_0','Phi_dense_0','Boundary_data_0','x_dense_half_0','Fb_ini_dense_half_0');

set(gca,'fontsize',20);
plot(x_dense_half_0,Fb_ini_dense_half_0,'s-.',x_dense_0,Phi_dense_0,'^-',x_dense_0,Boundary_data_0,'o-');
legend('v^3','f_\phi at infinity','f_\phi(v,x=0)','location','southwest');
xlabel('v');
xlim([min(x_dense_0) max(x_dense_0)]);
title(['given \phi (v>0) = v^3, u = 0']);
print(gcf,'-depsc2', 'data/Filtering2/u_convergence/N36.eps');
