clear;

% fb_ini_s.fun = inline('v.^2');
fb_ini_s.fun = inline('v');

u = 0;
ki = [];
ki = 4*[1:20];
ki = [80];
% ki = [ki,80];
% ki = [40,44,48,52,56,60,];

for k = 1:length(ki)
    N = ki(k);

%     [fb,X,W,H,fb_c] = boundary_coeff_damp(N,fb_ini_s);
    [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = boundary_coeff(N,fb_ini_s);
%     filename = 'data/damping/N_';
    filename = 'data/ToShow/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    save(filename,'eta','fb','g_tilde','phi','Boundary_data','X','W','H','fb_c');
%     save(filename,'fb','X','W','H','fb_c');
end