%% relative error
clear;

Ni = 4*[1:10]; Ni = [Ni,80]; Ni = Ni(:);
for k = 1:length(Ni)
    N = Ni(k);
    filename = 'data/damping/noFiltering_1/N_';
    filename = strcat(filename,num2str(N));
    filename = strcat(filename,'.mat');

    solution(k) = importdata(filename);
end

H_ref = solution(end).H;
W_ref = solution(end).W;

err = zeros(length(Ni)-1,1);
for k = 1:length(Ni) - 1;
    soln = H_ref(:,[1:Ni(k)-1,Ni(end):Ni(end)+Ni(k)-1])*solution(k).fb_c;
    error = soln - solution(end).fb;
    err(k) = error'*diag(W_ref)*error;
end

regression1 = polyfit(log(Ni(1:end-1)),log(err),[1]);
error_reg = polyval(regression1,log(Ni(1:end-1)));
error_reg = exp(error_reg);

h = figure(1);
set(gca,'fontsize',20);
% semilogy(Ni(1:end),error,'r.-.', Ni(1:end),error_reg(1:end),'bo-',...
%     Ni(1:end),error_optimal,'r.-.',Ni(1:end),error_reg_optimal(1:end),'bo-');
plot(log(Ni(1:end-1)),log(err(1:end)),'ro-',log(Ni(1:end-1)),log(error_reg));%,Ni,error_optimal,'r.-');
% legend('numerical error','regression','optimal error','regession','location','northeast');
% legend('numerical error','regression','optimal error','location','northeast');
legend('numerical error','regression','location','southwest');
xlabel('N');
ylabel('error');
% xlim([min(Ni)/1.2,2*max(Ni)]);
% ylim([1e-6 1e-2]);
title(['\phi_b (v>0) = 1']);
gtext(['slope = ',num2str(regression1(1),'%1.4f')]);
% text(20,10^(-3),['slope =', num2str(regression1(1))],'FontSize',16);
print(gcf,'-depsc2', 'data/damping/noFiltering_1/error_0_no-optimal.eps');
close(h)

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