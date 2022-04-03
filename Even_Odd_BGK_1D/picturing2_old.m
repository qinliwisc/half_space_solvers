%% relative error
clear;

% u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
u = 0.5;
ki = [2:5];
Ni = [4,8,12,16,20,24,28,32,34,36];
for k = 1:length(Ni)
    N = Ni(k);
    filename = 'data/u_convergence/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    
    solution(k) = importdata(filename);
end

for k = 2:length(Ni)
    ref = solution(k-1).Boundary_data;
    current = solution(k).Boundary_data;
    error(k) = norm(ref-current,2)/norm(ref,2);
end

set(gca,'fontsize',20);
loglog(Ni(2:end),error(2:end),'.-.');
% plot(log10(Ni(2:end)),log10(error(2:end)),'.-.');
% plot(ki,error,'.-.');
legend('error','location','northeast');
xlabel('N');
ylabel('error');
xlim([min(Ni)/1.2,2*max(Ni)]);ylim([min(error)/2,2*max(error)]);
title(['f_b (v>0) = v^3, u = 0.5']);
print(gcf,'-depsc2', 'data/u_convergence/error.eps');

%% absolute error
% set(gca,'fontsize',20);
% plot(log10(ki(3:end)),log10(error(3:end)),'.-.');
% legend('error','location','northwest');
% xlabel('log_{10}(N)');
% ylabel('error');
% title(['f_b (v>0) = \chi_0 + \chi_+, u = 0.5']);
% print(gcf,'-depsc2', 'data/u_convergence/error2.eps');

% clear;
% 
% % u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
% u = 0.5;
% ki = [2:5];
% % ki = [4,8,12,16,20,24,28,32,36];
% % ki = [4,8,12,16,20,22,24,26,28,30,32,34,36];
% for k = 2:length(ki)
%     N = 2^ki(k);
%     filename = 'data/u_convergence/N_';
%     filename = strcat(filename,num2str(N));
%     filename = strcat(filename,'.mat');
%     
%     solution(k) = importdata(filename);
%     error(k) = norm(solution(k).Boundary_data - solution(k).Fb_ini_dense);
% end
% 
% set(gca,'fontsize',20);
% plot(log10(ki),log10(error),'.-.');
% % plot(ki,error,'.-.');
% legend('error','location','northwest');
% xlabel('log_{10}(N)');
% ylabel('log_{10}(error)');
% title(['f_b (v>0) = \chi_0 + \chi_+, u = 0.5']);
% print(gcf,'-depsc2', 'data/u_convergence/error.eps');
% % 
% % set(gca,'fontsize',20);
% % plot(log10(ki(3:end)),log10(error(3:end)),'.-.');
% % legend('error','location','northwest');
% % xlabel('log_{10}(N)');
% % ylabel('error');
% % title(['f_b (v>0) = \chi_0 + \chi_+, u = 0.5']);
% % print(gcf,'-depsc2', 'data/u_convergence/error2.eps');