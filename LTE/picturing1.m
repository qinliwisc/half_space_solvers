clear

soln = load('data/index_threshold/N_60.mat');
mu = soln.X(soln.X>0);
[val, error] = Hfun(mu, 1, 0.0001);
real_soln = val/sqrt(3) - mu;
real_soln_neg = real_soln(end:-1:1);
real_soln = [mu;real_soln_neg];

figure(1)
set(gca,'fontsize',20);
plot(soln.X,soln.Boundary_data,'o',soln.X,real_soln);
legend('numerical solution','exact solution','location','southwest');
xlabel('v');
ylabel('f(x=0,v)');
ylim([-0.1 1]);
title(['\phi_b (v>0) = v']);
%print(gcf,'-depsc2', 'data/cos_threshold/compare_H.eps');