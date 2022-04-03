function [C,g_c] = ConstructC_xysplit(Nx,Ny,Nx_po, Ny_po,u,dampcoef)

%% discretization on v
[Xp,Wp] = half_hermquad_shift_symb(Nx-1,0);
Wp = Wp / sqrt(pi);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:Nx-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp - u; Xm = Xm - u; 

X = [Xp;Xm]; Wx = [Wp;Wm];
Hx = [Hp(:,1:Nx_po-1),Hp(:,1:Nx_po);...
    Hm(:,1:Nx_po-1),-Hm(:,1:Nx_po)];

vu = X + u;

[Y,Wy] = hermquad_symb(Ny-1);
[Y,ind] = sort(Y,'descend'); Wy = Wy(ind);
[Hy] = hermite_poly(Y, Ny-1) * (pi)^(1/4);
Wy = Wy/(Hy(:,1)'*diag(Wy)*Hy(:,1));

[yy,xx] = meshgrid(Y,vu);

%% polynomials basis
H_2d = kron(Hy(:,1:Ny_po),Hx);
W_2d = kron(Wy,Wx);
x_2d = kron(ones(length(Y),1),X);

chi(1).fun = inline('(vx.^2+vy.^2-2)/sqrt(2*sqrt(pi))','vx','vy');
chi(2).fun = inline('sqrt(2)*vy/sqrt(sqrt(pi))','vx','vy');
chi(3).fun = inline('(vx.^2+vy.^2+2*vx)/2/sqrt(sqrt(pi))','vx','vy');
chi(4).fun = inline('(vx.^2+vy.^2-2*vx)/2/sqrt(sqrt(pi))','vx','vy');

%% initial data
index_p = []; index_0 = [];

if (abs(u)<eps)
    index_p = [3];
    index_0 = [1,2];
elseif (u>0) 
    if (abs(u-sqrt(3/2))<eps)
        index_0 = [4];
        index_p = [1,2,3];
    elseif (u<sqrt(3/2))
        index_p = [1,2,3];
    elseif (u>sqrt(3/2))
        index_p = [1,2,3,4];
    end
elseif (u<0)
    if (abs(u+sqrt(3/2))<eps)
        index_0 = [3];
    elseif (u>-sqrt(3/2))
        index_p = [3];
    end
end

%% constructing matrix C
sizeC = length(index_p)+length(index_0);
C = zeros(sizeC,sizeC);
chip = []; chi0 = []; g0_c = []; gp_c = [];
for k=1:length(index_0)
    chi0_s = chi(index_0(k));
    chi0_d = feval(chi0_s.fun,xx,yy);
    chi0_d = reshape(chi0_d,length(X)*length(Y),1);
    chi0_d = x_2d.*chi0_d;
    chi0 = horzcat(chi0,chi0_d);
    [g0_c_d]=BoundaryData_2d_xysplit...
        (index_0(k), Nx, Ny, Nx_po, Ny_po, u, dampcoef);
    g0_c_d = g0_c_d(:);
    g0_c = horzcat(g0_c,g0_c_d);
%     g0_d = H_2d*g0_c_d;
%     g0 = horzcat(g0,g0_d);
end
for k=1:length(index_p)
    chip_s = chi(index_p(k));
    chip_d = feval(chip_s.fun,xx,yy);
    chip_d = reshape(chip_d,length(X)*length(Y),1);
    chip_d = x_2d.*chip_d;
    chip = horzcat(chip,chip_d);
    [gp_c_d]=BoundaryData_2d_xysplit...
        (index_p(k), Nx, Ny, Nx_po, Ny_po, u, dampcoef);
    gp_c_d = gp_c_d(:);
    gp_c = horzcat(gp_c,gp_c_d);
%     gp_d = H_2d*gp_c_d;
%     gp = horzcat(gp,gp_d);
end

chi = [chip,chi0];
g_c = [gp_c,g0_c];

C =chi'*diag(W_2d)*(H_2d*g_c);
return
