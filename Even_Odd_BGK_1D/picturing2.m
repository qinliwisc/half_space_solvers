%% relative error
clear;

% u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
u = 0;
Ni = 4*[1:6]; Ni = [Ni,36];
% Ni = [Ni,40];
Ni = Ni(:);%[4,8,12,16,20,24,28,32,34,36];
for k = 1:length(Ni)
    N = Ni(k);
    filename = 'data/withoutFiltering/u_convergence/u0/v_cube/old2/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');
    
    solution(k) = load(filename);
end

coeff = zeros(max(Ni)*2-1, length(Ni));
eta = zeros(2,length(Ni));
sigma = zeros(max(Ni)*2-1, length(Ni));
eta_c = 0.2; alpha = 2; p = 2;
for k=1:length(Ni)
    eta(:,k) = solution(k).eta;
    coeff(1:Ni(k)-1,k) = solution(k).coeff(1:Ni(k)-1);
    coeff(max(Ni):max(Ni)+Ni(k)-1,k) = solution(k).coeff(Ni(k):2* ...
                                                      Ni(k)-1);
end

ref = coeff(:, end);
ref_eta = eta(:,end);
error = zeros(length(Ni),1);
error_eta = zeros(length(Ni),1);
for k = 1:length(Ni)-1
    current = coeff(:, k);
    error(k) = norm(ref-current,2)/norm(ref,2);
end
regression1 = polyfit(log10(Ni(1:end-1)),log10(error(1:end-1)),[1]);
error_reg = polyval(regression1,log10(Ni(1:end-1)));
error_reg = 10.^(error_reg);
for k = 1:length(Ni)-1
    current_eta = eta(:,k);
    error_eta(k) = norm(ref_eta-current_eta,2)/norm(ref_eta,2);
end
regression2 = polyfit(log10(Ni(1:end-1)),log10(error_eta(1:end-1)),[1]);
error_eta_reg = polyval(regression2,log10(Ni(1:end-1)));
error_eta_reg = 10.^(error_eta_reg);

figure(1)
set(gca,'fontsize',20);
% loglog(Ni(1:end-1),error(1:end-1),'.-.');
plot(log10(Ni(1:end-1)),log10(error(1:end-1)),'ro-',log10(Ni(1:end-1)),log10(error_reg));
%semilogy(Ni(1:end-1),error_reg,'g.-.',Ni(1:end-1),error(1:end-1),'bo-');
legend('error at x= 0','regression','location','northeast');
xlabel('log(N)');
ylabel('log(error)');
title(['\phi_b (v>0) = v^3, u = 0']);
gtext(['slope = ',num2str(regression1(1))],'FontSize',16);
print(gcf,'-depsc2', 'data/withoutFiltering/u_convergence/u0/error_boundary.eps');

figure(2)
set(gca,'fontsize',20);
%  loglog(Ni(1:end-1),error_eta(1:end-1),'.-.');
plot(log10(Ni(1:end-1)),log10(error_eta(1:end-1)),'ro-',log10(Ni(1:end-1)),log10(error_eta_reg));
legend('error at x = \infty','regression','location','northeast');
xlabel('log(N)');
ylabel('log(error)');
title(['\phi_b (v>0) = v^3, u = 0']);
gtext(['slope =', num2str(regression2(1))],'FontSize',16);
print(gcf,'-depsc2', 'data/withoutFiltering/u_convergence/u0/error_infty.eps');

% %% relative error
% clear;
% 
% % u = [2,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
% u = 0;
% Ni = 2*[1:15]+2; Ni = [Ni,36];%[4,8,12,16,20,24,28,32,34,36];
% for k = 1:length(Ni)
%     N = Ni(k);
%     filename = 'data/u_convergence/u0/N_';
%     filename = strcat(filename,num2str(N));
%     filename = strcat(filename,'.mat');
%     
%     solution(k) = importdata(filename);
% end
% 
% 
% filenamex = 'prepare/x_'; filenamex = strcat(filenamex,num2str(N));
% filenamex = strcat(filenamex,'.dat');
% filenamew = 'prepare/w_'; filenamew = strcat(filenamew,num2str(N));
% filenamew = strcat(filenamew,'.dat');
% [Xp] = load(filenamex); Wp = load(filenamew);
% Wp = Wp / sqrt(pi);
% [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind);
% [Hp] = half_hermite_shift_poly(Xp, [0:1:2], 0) *(pi)^(1/4)/ sqrt(2);
% 
% Xm = -Xp; Wm = Wp; 
% [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
% [Hm] = half_hermite_shift_poly(abs(Xm), [0:1:2], 0) *(pi)^(1/4)/ sqrt(2);
% 
% chi(1).fun = inline('(2*v.^2-3)/sqrt(6)');
% chi0_p_data = feval(chi(1).fun,Xp);chi0_m_data = feval(chi(1).fun,Xm);
% chi(2).fun = inline('(sqrt(6)*v + 2*v.^2)/sqrt(6)');
% chip_p_data = feval(chi(2).fun,Xp);chip_m_data = feval(chi(2).fun,Xm);
% 
% for k=1:3
%     Ce_max(:,2) = Hp'*diag(Wp)*chi0_p_data + Hm'*diag(Wm)*chi0_m_data;
%     Ce_max(:,1) = Hp'*diag(Wp)*chip_p_data + Hm'*diag(Wm)*chip_m_data;
%     Co_max(:,2) = Hp'*diag(Wp)*chi0_p_data - Hm'*diag(Wm)*chi0_m_data;
%     Co_max(:,1) = Hp'*diag(Wp)*chip_p_data - Hm'*diag(Wm)*chip_m_data;
% end
% 
% coeff = zeros(127, length(Ni));
% for k=1:length(Ni)
%     len = length(solution(k).fb_c); 
%     coeff(1:(len-1)/2, k) = solution(k).fb_c(1:(len-1)/2);
%     coeff(63 + (1:(len-1)/2+1), k) = solution(k).fb_c((len-1)/2+1:end);
%     coeff(1:(len-1)/2, k) = coeff(1:(len-1)/2, k) - ...
%         solution(k).g_coeff(1:(len-1)/2);
%     coeff(63 + (1:(len-1)/2+1), k) = coeff(63 + (1:(len-1)/2+1), k) - ...
%         solution(k).g_coeff((len-1)/2+1:end);
%     coeff(1:3, k) = coeff(1:3, k) + Ce_max*solution(k).eta;
%     coeff(63+(1:3), k) = coeff(63+(1:3), k) + Co_max*solution(k).eta;
% end
% 
% ref = coeff(:, end);
% 
% for k = 1:length(Ni)-1
%     current = coeff(:, k);
% %     ref = solution(k-1).eta;
%     error(k) = norm(ref-current,2)/norm(ref,2);
% end
% 
% set(gca,'fontsize',20);
% % loglog(Ni(1:end-1),error(1:end-1),'.-.');
% semilogy(Ni(1:end-1),error(1:end),'.-.');
% legend('error','location','northeast');
% xlabel('N');
% ylabel('error');
% % xlim([min(Ni)/1.2,2*max(Ni)]);
% ylim([min(error)/1.2,2*max(error)]);
% title(['f_b (v>0) = v^3, u = 0']);
% print(gcf,'-depsc2', 'data/u_convergence/u0/error_boundary.eps');

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