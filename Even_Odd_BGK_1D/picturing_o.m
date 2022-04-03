clear;

% fb_ini_s.fun = inline('(2*v.^2-3)/sqrt(6)');

fb_ini_s.fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');

% fb_ini_s.fun = inline('(sqrt(6)*v - 2*(v).^2)/sqrt(6)');

N = 16;
u = [-sqrt(1.5),-0.5,0,0.5,sqrt(1.5),2];
for k = 1:length(u)
%     [x_dense,Fb_ini_dense,Phi_dense,vu,Boundary_data]=boundary(N,u(k),fb_ini_s);
        [x_dense,Fb_ini_dense,Phi_dense,Boundary_data_dense]=boundary(N,u(k),fb_ini_s);
%     [Fb_ini,Phi,Boundary_data,vu,X,Xp]=boundary(N,u(k),fb_ini_s);
%     [x_dense,Fb_ini,Fb]=boundary(N,u(k),fb_ini_s);
    x_dense_half = x_dense(x_dense>-u(k));
    Fb_ini_dense_half = Fb_ini_dense(x_dense>-u(k));
%     plot(x_dense_half,Fb_ini_half,'.-.',x_dense,Fb,'o-');
%     legend('real: 0 mode','recover','location','southeast');
%     xlabel('v');
%     ylabel('f_b');
    plot(x_dense_half,Fb_ini_dense_half,'.-.',x_dense,Phi_dense,'^-',x_dense,Boundary_data_dense,'o-');
%     plot(Xp,Fb_ini,'.-.',X,Phi,'r',X,Boundary_data,'o-');
    legend('real: \chi_+ mode','f at infinity','recover, full','location','northwest');
    xlabel('v');
    ylabel('f(v)');
    title(['given f_b (v>0) = \chi_+ mode, u = ',num2str(u(k))]);
    filename = 'data/x-plus/u=';
    filename = strcat(filename,num2str(u(k)));
    filename = strcat(filename,'.eps');
    print('-depsc', filename);
end