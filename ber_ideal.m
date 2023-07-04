%% ber_ideal

%% parameters
SNRdB = [0:1:30];
Nu = par.Nu;
Nt = par.Nt;
Bw = par.Bw;
N_carrier = 128;
fs = par.fs;
dt = par.dt;
N_base = 1;
T_diff = round(dt*Bw);
% Bw = 20e6;
% N_carrier = 64;
N_ofdm = 10;
N_cp = N_carrier/4;
T = (N_carrier+N_cp)*N_ofdm;
Nloop = 200;
Ms = [16];
Ber_ideal = zeros(length(SNRdB),length(Ms));
sigPow = 0;
pwcons = 0;
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
         %% Ideal 
         Ber_ideal(n) =  Ber_ideal(n) + IdeaBer(data, s_ofdm, M, Nu, T, sigma, N_carrier, N_cp);
      end
      Ber_ideal(n) = Ber_ideal(n)/Nloop;
   end
   Ber_ideal2(:,l) = Ber_ideal;
end

figure;
semilogy(SNRdB,Ber_ideal2(:),'-o','linewidth',2), grid on
% 
legend('Ideal(AWGN)','Location','Southwest','fontsize',16)
xlim([1 30])
%%
save Ideal_ber.mat Ber_ideal2;