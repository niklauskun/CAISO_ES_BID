% x= [0	0.1	1	10	20	100];
% y1= [18757	18755	18740	18529	18066	13237];
% y2= [14948	14948	14987	15490	15490	8575];
% plot(x,y1,'LineWidth',5) 
% hold on 
% plot(x,y2,'LineWidth',5)
% hold off
% legend('Linear ES', 'Nonlinear ES')
% xlabel('\sigma')
% ylabel('Profit ($)')
% ylim([0 20000])
% set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

% x= [1 5 10 20];
% y1= [16908	18757	18666	17847];
% y2= [10042	14948	14740	14821];
% plot(x,y1,'LineWidth',5) 
% hold on 
% plot(x,y2,'LineWidth',5)
% hold off
% legend('Linear ES', 'Nonlinear ES')
% xlabel('Number of segments')
% ylabel('Profit ($)')
% ylim([0 20000])
% set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

cbid = csvread('./Data/WALNUT_2016_Charge.csv');
a = cbid(4,:) - cbid(5,:)>0;


