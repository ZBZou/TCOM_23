%%BER plot

%%
c = load('BER_space_Time_4_6.mat');
c_ideal = load('Ideal_ber.mat');
c_dpc = load('DPC_ber.mat');
SNRdB = c.SNRdB;
Ber_hogmt =  c.Ber2;
ideal_Ber = c_ideal.Ber_ideal2;
DPC_ber = c_dpc.BER';
%%
figure;
semilogy(SNRdB, DPC_ber(1:length(SNRdB),1),'-o',...
   SNRdB,Ber_hogmt(:,1),'-o',SNRdB,Ber_hogmt(:,2),'-o',...
   SNRdB,Ber_hogmt(:,3),'-o',SNRdB,Ber_hogmt(:,4),'-o',...
   SNRdB,Ber_hogmt(:,5),'-o',SNRdB,ideal_Ber(1:length(SNRdB),1),'-o','linewidth',2), grid on
% 
legend('DPC','HOGMT-Precoding: 90%','HOGMT-Precoding: 95%','HOGMT-Precoding: 98%',...
   'HOGMT-Precoding: 99%','HOGMT-Precoding: 100%',...
   'Ideal(AWGN)','Location','Southwest','fontsize',16)
xlim([1 20])
xlabel ('SNR (dB)')
ylabel('BER (dB)')
set(gca, 'fontsize', 20)
[h, wd, ht] = tightfig();
filename = 'Ber_v120';
save_hst = 0;
if save_hst == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end