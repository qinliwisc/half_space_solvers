% function hermite_2d
clear;

N = 4; u = 0.5;

[Xp,Wp] = half_hermquad_shift_symb(N-1,0);
[Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
[Hp] = half_hermite_shift_poly(Xp, [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xm = -Xp; Wm = Wp; 
[Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
[Hm] = half_hermite_shift_poly(abs(Xm), [0:1:N-1], 0) *(pi)^(1/4)/ sqrt(2);

Xp = Xp-u;
Xm = Xm-u; 


X = [Xp;Xm]; W = [Wp;Wm]; W = W / sqrt(pi);

vu = X + u;


%% 
H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

H_2d = kron(H,H);
W_2d = kron(W,W);
x_2d = kron(X,ones(size(X)));
y_2d = kron(ones(size(X)),X);

% return