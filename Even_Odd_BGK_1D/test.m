clear;

fb_ini_s_0.fun = inline('(2*v.^2-3)/sqrt(6)');

fb_ini_s_p.fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');

fb_ini_s_m.fun = inline('(sqrt(6)*v - 2*(v).^2)/sqrt(6)');

N = 16;
u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];

for k=3:-1:-2
    j = -k+4;
    [x_dense_0,Fb_ini_dense_0,Phi_dense_0,Boundary_data_0]=boundary(N,u(j),fb_ini_s_0,8);
    [x_dense_p,Fb_ini_dense_p,Phi_dense_p,Boundary_data_p]=boundary(N,u(j),fb_ini_s_p,8);
    [x_dense_m,Fb_ini_dense_m,Phi_dense_m,Boundary_data_m]=boundary(N,u(j),fb_ini_s_m,8);

    x_dense_half_0 = x_dense_0(x_dense_0>-u(j));
    Fb_ini_dense_half_0 = Fb_ini_dense_0(x_dense_0>-u(j));
    x_dense_half_p = x_dense_p(x_dense_p>-u(j));
    Fb_ini_dense_half_p = Fb_ini_dense_p(x_dense_p>-u(j));
    x_dense_half_m = x_dense_m(x_dense_m>-u(j));
    Fb_ini_dense_half_m = Fb_ini_dense_m(x_dense_m>-u(j));

    filename = 'data/Filtering2/u';
    filename = strcat(filename,num2str(k));
    filename_0 = strcat(filename,'/zero.mat');
    filename_p = strcat(filename,'/plus.mat');
    filename_m = strcat(filename,'/minus.mat');
    save(filename_0,'x_dense_0','Fb_ini_dense_0','Phi_dense_0','Boundary_data_0','x_dense_half_0','Fb_ini_dense_half_0');
    save(filename_p,'x_dense_p','Fb_ini_dense_p','Phi_dense_p','Boundary_data_p','x_dense_half_p','Fb_ini_dense_half_p');
    save(filename_m,'x_dense_m','Fb_ini_dense_m','Phi_dense_m','Boundary_data_m','x_dense_half_m','Fb_ini_dense_half_m');

end