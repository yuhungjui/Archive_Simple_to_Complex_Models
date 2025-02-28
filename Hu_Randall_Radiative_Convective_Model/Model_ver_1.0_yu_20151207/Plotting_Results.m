clear; close all; clc;
% =========================================================================
% 
% Radiative-Convective Systems Idealized Model.
% Formulations based on Hu and Randall 1995: Low-Frequency Oscillations in
% Radiative-COnvective Systems. Part II: An Idealized Model.
% 
% Plotting Results.
% 
% =========================================================================
%% Load Results:

load ./Output/Test_1.mat

% =========================================================================
%% Plotting:

figure(1)

subplot(3,2,1)
plot(tday(1:smax),SM(1:smax),'LineWidth',1.2,'Color','r')

% axis([0 max(tday) -inf inf])

subplot(3,2,2)
plot(tday(1:smax),SBt(1:smax),'LineWidth',1.2,'Color','r')

% axis([0 max(tday) -inf inf])

subplot(3,2,3)
plot(tday(1:smax),qM(1:smax),'LineWidth',1.2,'Color','b')

% axis([0 max(tday) -inf inf])

subplot(3,2,4)
plot(tday(1:smax),qBt(1:smax),'LineWidth',1.2,'Color','b')

% axis([0 max(tday) -inf inf])

subplot(3,2,5)
plot(tday(1:smax),Mc(1:smax),'LineWidth',1.2,'Color','k')

% axis([0 max(tday) -inf inf])

subplot(3,2,6)
plot(tday(1:smax),PR(1:smax),'LineWidth',1.2,'Color','g')

% axis([0 max(tday) -inf inf])

% a1 = gca;
% set(a1,'Box','on');
% set(a1,'TickDir','out');
% set(a1,'Linewidth',1.2);
% set(a1,'FontName','FixedWidth');
% set(a1,'FontSize',14);
% set(a1,'FontWeight','Bold');
% set(a1,'XScale','linear');
% % set(a1,'XTick',[]);
% % set(a1,'XTickLabel',{});
% set(a1,'XMinorTick','off');
% set(a1,'XMinorGrid','off');
% % set(a1,'YTick',[]);
% % set(a1,'YTickLabel',{});
% set(a1,'YMinorTick','off');
% set(a1,'YMinorGrid','off');
% % a1.YDir = 'Reverse';
% 
% xlabel('\bf{Day}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% % ylabel('\bf{}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
% 
% % tit2 = title('','FontName','Helvetica','FontSize',14,'FontWeight','bold');
% 
% % leg = legend('OND (10/01-12/31)','MJO1 (10/15-10/31)',2);
% % set(leg,'FontName','Helvetica')
% % set(leg,'FontSize',14)
% % set(leg,'FontWeight','normal')
% % set(leg,'Location','northwest')

% =========================================================================
%% Save Figure:

% set(gf1,'Color',[1,1,1]);
% figname = ['./Results'];
% print(gf1,'-dpng','-r600',figname);

% =========================================================================
