clear;

fb_ini_s_0.fun = inline('(2*v.^2-3)/sqrt(6)');

fb_ini_s_p.fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');

fb_ini_s_m.fun = inline('(sqrt(6)*v - 2*(v).^2)/sqrt(6)');

% fb_ini_s.fun = inline('0.4*(2*v.^2-3)/sqrt(6)+0.6*(sqrt(6)*v + 2*v.^2)/sqrt(6)');
fb_ini_s.fun = inline('v.^3');

% u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
u = 0;
ki = 40;
% ki = 4*[1:9];
% ki = [40,44,48,52,56,60,];

for k = 1:length(ki)
    N = ki(k);
    [fb_c,g_coeff,eta,phi_coeff,coeff,X,Boundary_data] = boundary_coeff(N,u,fb_ini_s);
%     [fb_c,g_coeff,eta,index_total] = boundary_coeff(N,u,fb_ini_s);
%     
    filename = 'data/u_convergence/u0/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    save(filename,'fb_c','g_coeff','eta','phi_coeff','coeff','X','Boundary_data');
end