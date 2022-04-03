% function legendre_2d
clear;

N = 4;

% 1D
[Xp, Wp] = half_legendre_quad(N-1);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
Hp = half_legendre_poly(Xp,[0:1:N-1]);

Xm = -Xp; Wm = Wp;
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
X = [Xp;Xm]; W = [Wp;Wm]; W = W/2;

H_2d = kron(H,H);
W_2d = kron(W,W);
x_2d = kron(X,ones(size(X)));
y_2d = kron(ones(size(X)),X);


% return