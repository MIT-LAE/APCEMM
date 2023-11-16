clear all;
close all;

T_list = 217-10:1:217+20;
RH_list = 0:10:140;

T_evals_matrix = readmatrix('APCEMM-sweep-evaluations-T.csv')
T_evals = T_evals_matrix(2:end,end);

RH_evals_matrix = readmatrix('APCEMM-sweep-evaluations-RH.csv')
RH_evals = RH_evals_matrix(2:end,end);

Times_matrix = readmatrix('APCEMM-sweep-times.csv')
Times = Times_matrix(2:end,end);




figure(1)
% set figure position and size:
set(gcf,'position',[160 200 500 400]) % 350

%semilogy(T_list, T_evals)
plot(T_list, T_evals)




figure(2)
% set figure position and size:
set(gcf,'position',[160 200 500 400]) % 350

%semilogy(RH_list, RH_evals)
plot(RH_list, RH_evals)




% hold on
% scatter(x, y, 100, 'k',Marker='x',MarkerFaceColor='k')
% hold off
% 
% box off
% 
% xlim([0,4.5]);
% % xticks(linspace(210,260,3));
% xlabel('Ease of Implementation Rank');
% 
% ylim([0,4.5]);
% % yticks(0:50:150);
% ylabel('Effectiveness Rank')
% 
% %set(gca,'XTickLength',[0 0])
% 
% % keep position and size when printing:
% set(gcf,'PaperPositionMode','auto')
% 
% % set fonts and frame:
% set(gca,'Fontn','Arial','FontSize',16,'linewidth',1)
% 
% % For a vectorial figure (for latex):
% print -deps2c 01-Aviation-Emissions.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 01-Aviation-Emissions.png