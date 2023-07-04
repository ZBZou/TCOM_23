%% BER
clear all;
close all;
%% load
c = load ('h_d10_v120.mat');
par = load('par_d10_v120.mat');
%% parameters
Nu = par.Nu;
Nt = par.Nt;
Bw = par.Bw;
N_carrier = 64;
fs = par.fs;
dt = par.dt;
N_base = 1;
T_diff = round(dt*Bw);
% Bw = 20e6;
% N_carrier = 64;
N_ofdm = 10;
N_cp = N_carrier/4;
Nu = size(c.h,1);
T = (N_carrier+N_cp)*N_ofdm;

Nloop = 5;
ep = [0.9, 0.95, 0.98, 0.99, 1];
Nep = length(ep);
SNRdB = [0:1:20];
Ber = zeros(length(SNRdB),Nep);
Ber_dpc = zeros(length(SNRdB),1);
Ber_ideal = zeros(length(SNRdB),1);

Ms = [16];
P = zeros(T,Nep);
Ber2 = zeros(length(SNRdB),Nep, length(Ms));
Ber_dpc2 = zeros(length(SNRdB),1, length(Ms));
Ber_ideal2 = zeros(length(SNRdB),length(Ms));
sigPow = 0;
pwcons = 0;
pw = 1;

%% interpolate
hc = zeros(Nu,Nt,size(c.h,3),T);
for u1 = 1:Nu
   for u2 = 1:Nt
      for t1 = 1:size(c.h,3)
         for t2 = 1:T
            t_temp = round(t2/T_diff)+1;
            hc(u1,u2,t1,t2) = c.h(u1,u2,t1,t_temp);
         end
      end
   end
end

%% CSI 
% Hc = c.h(1:Nu,1:Nu,:,1:T);

Hc = permute(hc, [1,4,2,3]);

Hst = zeros(Nu,T,Nu,T);

for u1 = 1:Nu
    for u2 = 1:Nu
        for t1 = 1:T
            for t2 = t1:t1+size(Hc,4)-1
                if t2 > T
                    continue;
                end
                Hst(u1,t1,u2,t2) = Hc(u1,t1,u2,t2-t1+1);
            end
        end
    end
end

%% Loop

for l = 1:length(Ms)
   M = Ms(l);
   K = log2(M);
   %% data
   data = randi([0 1],Nu,N_carrier*N_ofdm*K);
   s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
   s = s';
   %% OFDM
   s_temp = reshape(s, Nu, N_carrier, []);
   s_temp = permute(s_temp, [2,3,1]);
   s_ofdm = ofdmmod(s_temp, N_carrier, N_cp);
   s_ofdm = permute(s_ofdm, [2,1]);
   
   %% OFDM_HOGMT
   s_ofdm_hogmt = ofdmmod(s_temp, N_carrier, 0);
   s_ofdm_hogmt = permute(s_ofdm_hogmt, [2,1]);
   s_ofdm_hogmt = [s_ofdm_hogmt,s_ofdm_hogmt(:,1:1/4*size(s_ofdm_hogmt,2))];
   
   for n = 1:length(SNRdB)
      sigma2 = 0.5*10^(-SNRdB(n)/10); sigma = sqrt(sigma2);
      sq2 = square(2);
      for m = 1:Nloop
         display(l,['M-QAM']);
         display(n,['SNR']);
         display(m,['Loop']);
         
         %% HOGMT

         X_hogmt = HOGMT(s_ofdm_hogmt, Hst, ep);
         %% power constriant
%          if pwcons == 1
%             p_hogmt = sum(sum(abs(X_hogmt).^2));
%             X_hogmt = X_hogmt * sqrt(pw/p_hogmt);
%          end
         %% receive Rx

         Rx = zeros(Nu, T, Nep);
         for i = 1:Nu
            for j = 1:T
                for ep_i = 1:Nep
                   Rx(i,j,ep_i) = sum(sum(squeeze(Hst(i,j,:,:)) .* X_hogmt(:,:,ep_i)));
                end
            end
         end 
         %% AWGN       
         for i = 1:Nu
            for j = 1:T
                for ep_i = 1:Nep
                   sigPow = sum(sum(abs(Rx(:,:,ep_i)),1),2)/(size(Rx(:,:,ep_i),1)*size(Rx(:,:,ep_i),2));
                   Rx(i,j,ep_i) = Rx(i,j,ep_i) + sigPow*sigma*(randn()+randn()*1i)/sqrt(2);
                end
            end
         end 
         %% Normalization
%          if pwcons == 1
%             Rx = Rx * sqrt(p_hogmt/pw);
%          end
         %% demod_HOGMT
         for ep_i = 1:Nep
             temp_s = Rx(:,1:4/5*size(Rx,2),ep_i);
             r_temp = permute(temp_s, [2,1]);
             s_demod = ofdmdemod(r_temp, N_carrier, 0);
             s_hat = reshape(permute(s_demod, [3,1,2]),Nu,[]);
  
             
             %% ber
            
             Ber(n,ep_i) = Ber(n,ep_i) + myber(data,s_hat,M);
         end
         %% Ideal 
         Ber_ideal(n) =  Ber_ideal(n) + IdeaBer(data, s_ofdm, M, Nu, T, sigma, N_carrier, N_cp);
         %% DPC
%          ne_dpc = ne_DPC(Nu, T, s, data, Hst, sigma, M, taps_max,mode);
% %          ne_dpc = 0;
%          Ber_dpc(n) = Ber_dpc(n) + ne_dpc/lengh(data);
      end
      Ber(n,:) = Ber(n,:)/Nloop;
%       Ber_dpc(n,:) = Ber_dpc(n,:)/Nloop;
      Ber_ideal(n) = Ber_ideal(n)/Nloop;
   end
   Ber2(:,:,l) = Ber;
%    Ber_dpc2(:,:,l) = Ber_dpc;
   Ber_ideal2(:,l) = Ber_ideal;
end


%% save
save_data = 0;
if save_data == 1
   save BER_hogmt_d10_v120_2.mat Ber2 Ber_ideal2 ep Ms SNRdB;
end

%%
figure;
semilogy(SNRdB,Ber2(:,1),'-o',SNRdB,Ber2(:,2),'-o',...
   SNRdB,Ber2(:,3),'-o',SNRdB,Ber2(:,4),'-o',...
   SNRdB,Ber2(:,5),'-o',SNRdB,Ber_ideal2(:,1),'-o','linewidth',2), grid on
% 
legend('HOGMT-Precoding: 90%','HOGMT-Precoding: 95%','HOGMT-Precoding: 98%',...
   'HOGMT-Precoding: 99%','HOGMT-Precoding: 100%',...
   'Ideal(AWGN)','Location','Southwest','fontsize',16)
xlim([1 30])
% semilogy(SNRdB,Ber_dpc2(:,1),'-o',SNRdB,Ber2(:,1),'-o',SNRdB,Ber2(:,2),'-o',...
%    SNRdB,Ber2(:,3),'-o',SNRdB,Ber_ideal2(:,1),'-o','linewidth',2), grid on
% 
% legend('DPC','HOGMT-Precoding: \epsilon {=} 1','HOGMT-Precoding: \epsilon {=} 10^{-1}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-2}','Ideal(AWGN)','Location','Southwest','fontsize',16)

% semilogy(SNRdB,Ber_dpc2(:,1),'-o',SNRdB,Ber2(:,1),'-o',SNRdB,Ber2(:,2),'-o',...
%    SNRdB,Ber2(:,3),'-o',SNRdB,Ber_ideal2(:,1),'-o','linewidth',2), grid on
% 
% legend('DPC with equalization','HOGMT-Precoding: \epsilon {=} 10^{-2}','HOGMT-Precoding: \epsilon {=} 10^{-3}',...
%    'HOGMT-Precoding: \epsilon {=} 10^{-4}','Ideal(AWGN)','Location','Southwest','fontsize',16)


% semilogy(SNRdB,Ber2(:,1,1),'-o',SNRdB,Ber2(:,1,2),'-o',SNRdB,Ber2(:,1,3),'-o',...
%    SNRdB,Ber2(:,1,4),'-o','linewidth',2), grid on
% 
% 
% legend('HOGMT-Precoding:BPSK','HOGMT-Precoding:QPSK',...
%    'HOGMT-Precoding:16QAM','HOGMT-Precoding:64QAM','Location','Southwest','fontsize',16)
xlabel ('SNR (dB)')
ylabel('BER (dB)')
set(gca, 'fontsize', 20)
[h, wd, ht] = tightfig();
filename = 'Ber_space_time';
save_hst = 0;
if save_hst == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end
% if mode == 0
%    print -opengl -dpdf -r600 Ber_space_time_temp.pdf
% elseif mode == 1
%    print -opengl -dpdf -r600 Ber_space_temp.pdf   
% end
%% plot H_s
% Hs = squeeze(Hst(:,1,:,1));
% mesh((abs(Hs)));
% xlabel('u^\prime');
% xlabh = get(gca,'XLabel');
% ylabel('u');
% ylabh = get(gca,'YLabel');
% set(gca,'fontsize',30)
% 
% [h, wd, ht] = tightfig();
% print -opengl -dpdf -r600 Hs.pdf
