clear;

N = 80;

fb_ini_s1.fun = inline('v');
fb_ini_s2.fun = inline('v');

sigma = cell(3,1);
alpha = cell(3,1);
sigma{1} = [0.5,0.5;0.5,0.5];
sigma{2} = [0.25,sqrt(0.25*0.75);sqrt(0.25*0.75),0.75];
sigma{3} = [0.25,0.25;0.75,0.75];
alpha{1} = [0.8,0.1,0.1];
alpha{2} = [1,0,0];
alpha{3} = [0.1,0.45,0.45];

for km = 1:3
    for kn = 1:3
        [eta,fb,g_tilde,phi,Boundary_data,X,W,H,fb_c] = ...
    boundary_coeff_Maxwell_multi(N,fb_ini_s1,fb_ini_s2,sigma{km},alpha{kn},(km-1)*3+kn);
    end
end