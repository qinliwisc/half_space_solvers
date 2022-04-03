function [C,g_c_d] = ConstructC(N,X,W,H,fb_ini_s1,fb_ini_s2)

%% discretization on v
vu = X;

%% constructing matrix C
sizeC = 1;%length(index_p)+length(index_0);
C = zeros(sizeC,sizeC);

g_c_d = BoundaryData_multi(N,1,fb_ini_s1,fb_ini_s2);
H_tilde = [H,zeros(size(H));zeros(size(H)),H];
g_d = H_tilde*g_c_d;

vuchi = vu.*feval(chi.fun,X); vuchi = vuchi(:); vuchi = [vuchi;vuchi];

C = vuchi'*diag([W;W])*g_d;
return