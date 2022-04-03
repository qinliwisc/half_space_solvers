function [C,g_c_d] = ConstructC(N,u,index_p,index_0,X,W,H,fb_ini_s)

%% discretization on v
vu = X + u;

chi(1).fun = inline('(2*v.^2-3)/sqrt(6)');
chi(2).fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
chi(3).fun = inline('(sqrt(6)*v - 2*v.^2)/sqrt(6)');

index = [index_p,index_0];

%% constructing matrix C
sizeC = length(index_p)+length(index_0);
C = zeros(sizeC,sizeC);
vuchip = []; vuchi0 = [];

g_c_d = BoundaryData_c(N,u,index,X,W,H,fb_ini_s);
g_d = H*g_c_d;
for k=1:length(index_0)
    chi0_s = chi(index_0(k));
    vuchi0_d = vu.*feval(chi0_s.fun,X); vuchi0_d = vuchi0_d(:);
    vuchi0 = horzcat(vuchi0,vuchi0_d);
end
for k=1:length(index_p)
    chip_s = chi(index_p(k));
    vuchip_d = vu.*feval(chip_s.fun,X); vuchip_d = vuchip_d(:);
    vuchip = horzcat(vuchip,vuchip_d);
end

vuchi = horzcat(vuchip,vuchi0);

C = vuchi'*diag(W)*g_d;
return