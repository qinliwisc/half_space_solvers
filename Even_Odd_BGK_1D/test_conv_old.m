clear;

fb_ini_s_0.fun = inline('(2*v.^2-3)/sqrt(6)');

fb_ini_s_p.fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');

fb_ini_s_m.fun = inline('(sqrt(6)*v - 2*(v).^2)/sqrt(6)');

% fb_ini_s.fun = inline('0.4*(2*v.^2-3)/sqrt(6)+0.6*(sqrt(6)*v + 2*v.^2)/sqrt(6)');
fb_ini_s.fun = inline('v.^3');

% u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
u = 0.5;
ki = [4,8,12,16,20,24,28,32,34,36];
for k = 1:length(ki)
    N = ki(k);
    [x_dense,Fb_ini_dense,Phi_dense,Boundary_data]=boundary(N,u,fb_ini_s);

    x_dense_half = x_dense(x_dense>-u);
    Fb_ini_dense_half = Fb_ini_dense(x_dense>-u);

    filename = 'data/u_convergence/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    save(filename,'x_dense','Fb_ini_dense','Phi_dense','Boundary_data','x_dense_half','Fb_ini_dense_half');
end