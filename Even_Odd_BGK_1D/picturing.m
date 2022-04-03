clear;

u = [2];%,sqrt(1.5),0.5,0,-0.5,-sqrt(1.5)];
for k = 3%:-1:-2
    filename = 'data/Filtering2/u';
    filename = strcat(filename,num2str(k));
    filename_0 = strcat(filename,'/zero.mat');
    filename_p = strcat(filename,'/plus.mat');
    filename_m = strcat(filename,'/minus.mat');

    zero = importdata(filename_0);
    plus = importdata(filename_p);
    minus = importdata(filename_m);

    filename_0 = strcat(filename,'/zero.eps');
    filename_p = strcat(filename,'/plus.eps');
    filename_m = strcat(filename,'/minus.eps');
    

    j = -k+4;j=1;
    set(gca,'fontsize',20);
%     [x_dense_0,Fb_ini_dense_0,Phi_dense_0,Boundary_data_0]=boundary(N,u(j),fb_ini_s_0);
    plot(zero.x_dense_half_0,zero.Fb_ini_dense_half_0,'s-.',zero.x_dense_0,zero.Phi_dense_0,'^-',zero.x_dense_0,zero.Boundary_data_0,'o-');
%     plot(zero.x_dense_half_0,zero.Fb_ini_dense_half_0,'.-.',zero.x_dense_0,zero.Phi_dense_0,'r',zero.vu_0,zero.Boundary_data_0,'o-');
    legend('\chi_0','f_\phi at infinity','f_\phi(v,x=0)','location','southwest');
    xlabel('v');
%     ylabel('f(v)');
    xlim([min(zero.x_dense_0) max(zero.x_dense_0)]);
    title(['given \phi (v+u>0) = \chi_0, u = ',num2str(u(j))]);
    print(gcf,'-depsc2', filename_0);
    
    set(gca,'fontsize',20);
    plot(plus.x_dense_half_p,plus.Fb_ini_dense_half_p,'s-.',plus.x_dense_p,plus.Phi_dense_p,'^-',plus.x_dense_p,plus.Boundary_data_p,'o-');
%     plot(plus.x_dense_half_p,plus.Fb_ini_dense_half_p,'.-.',plus.x_dense_p,plus.Phi_dense_p,'r',plus.vu_p,plus.Boundary_data_p,'o-');
    legend('\chi_+','f_\phi at infinity','f_\phi(v,x=0)','location','northwest');
    xlabel('v');
%     ylabel('f(v)');
    xlim([min(plus.x_dense_p) max(plus.x_dense_p)]);
    title(['given \phi (v+u>0) = \chi_+, u = ',num2str(u(j))]);
    print(gcf,'-depsc2', filename_p);
    
    
    set(gca,'fontsize',20);
    plot(minus.x_dense_half_m,minus.Fb_ini_dense_half_m,'s-.',minus.x_dense_m,minus.Phi_dense_m,'^-',minus.x_dense_m,minus.Boundary_data_m,'o-');
%     plot(minus.x_dense_half_m,minus.Fb_ini_dense_half_m,'.-.',minus.x_dense_m,minus.Phi_dense_m,'r',minus.vu_m,minus.Boundary_data_m,'o-');
    legend('\chi_-','f_\phi at infinity','f_\phi(v,x=0)','location','southwest');
    xlabel('v');
%     ylabel('f(v)');
    xlim([min(minus.x_dense_m) max(minus.x_dense_m)]);
    title(['given \phi (v+u>0) = \chi_-, u = ',num2str(u(j))]);
    print(gcf,'-depsc2', filename_m);
end
% 
% zero = importdata('data/u2/zero.mat');
% plus = importdata('data/u2/plus.mat');
% minus = importdata('data/u2/minus.mat');
% 
% figure(1)
% set(gca,'fontsize',20);
% suptitle(['given f_b (v>0) = \chi_- mode, u = 2']);
% subplot(3,1,1);
% plot(zero.x_dense_half_0,zero.Fb_ini_dense_half_0,'.-.',zero.x_dense_0,zero.Phi_dense_0,'r',zero.vu_0,zero.Boundary_data_0,'o-');
% legend('real: \chi_0 mode','f at infinity','recover, full','location','northwest');
% xlabel('v');
% ylabel('f(v)');
% subplot(3,1,2);
% plot(x_dense_half_p,Fb_ini_dense_half_p,'.-.',x_dense_p,Phi_dense_p,'r',vu_p,Boundary_data_p,'o-');
% legend('real: \chi_+ mode','f at infinity','recover, full','location','northwest');
% xlabel('v');
% ylabel('f(v)');
% subplot(3,1,3);
% plot(x_dense_half_m,Fb_ini_dense_half_m,'.-.',x_dense_m,Phi_dense_m,'r',vu_m,Boundary_data_m,'o-');
% legend('real: \chi_- mode','f at infinity','recover, full','location','northwest');
% xlabel('v');
% ylabel('f(v)');
% filename = 'data/u2/u2.eps';
% print(gcf,'-depsc2',filename);