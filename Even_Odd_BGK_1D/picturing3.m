%% relative error
clear;

% u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
u = 0;
Ni = 4*[1:9];% Ni = [Ni,36];%[4,8,12,16,20,24,28,32,34,36];
for k = 1:length(Ni)
    N = Ni(k);
    filename = 'data/u_convergence/u0/v_qua/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    
    solution(k) = importdata(filename);
end

ref = solution(end).eta;

for k = 2:length(Ni)
    current = solution(k).eta;
    ref = solution(k-1).eta;
    error(k) = norm(ref-current,2)/norm(ref,2);
end

set(gca,'fontsize',20);
% loglog(Ni(1:end-1),error(1:end-1),'.-.');
semilogy(Ni(2:end),error(2:end),'.-.');
legend('error','location','northeast');
xlabel('N');
ylabel('error');
% xlim([min(Ni)/1.2,2*max(Ni)]);
ylim([min(error)/1.2,2*max(error)]);
title(['f_b (v>0) = v^4, u = 0']);
print(gcf,'-depsc2', 'data/u_convergence/u0/v_qua/error_infty.eps');