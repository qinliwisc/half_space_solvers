function [C,g_c_d] = ConstructC(N,X,W,H,fb_ini_s)

%% discretization on v
vu = X;

%% constructing matrix C
sizeC = 1;%length(index_p)+length(index_0);
C = zeros(sizeC,sizeC);

g_c_d = BoundaryData_c(N,1,fb_ini_s);
g_d = H*g_c_d;

vuchi = vu.*feval(chi.fun,X); vuchi = vuchi(:);

C = vuchi'*diag(W)*g_d;
return