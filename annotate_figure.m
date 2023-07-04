%%figure
close all;

openfig('boxplot.fig')
% 
% legend('DPC with equalization','HOGMT-Precoding: \epsilon {=} 10^{-2}','HOGMT-Precoding: \epsilon {=} 10^{-3}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-4}','Ideal(AWGN)','Location','Southwest','fontsize',16)
% 
% grid on;
% box on;
set(gca, 'fontsize', 20)
% set(0,'defaultTextFontSize', 20)                      	% Default Font Size
% set(0,'defaultAxesFontSize', 20)  
% set(0,'defaultAxesFontSize',20)
% set(gca,'XTick',1:5:21);
% set(gca,'XTickLabel',(1:5:21));
% h = legend;
% set(h,'FontSize',20);
filename = 'boxplot';
[h, wd, ht] = tightfig();
%%
save_file = 1;
if save_file == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end
% print -opengl -dpdf -r600 Ber_multi_mod.pdf


% legend('DPC','HOGMT-Precoding: \epsilon {=} 10^{-1}','HOGMT-Precoding: \epsilon {=} 5{\times}10^{-2}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-2}','Ideal(AWGN)','Location','Southwest','fontsize',16)