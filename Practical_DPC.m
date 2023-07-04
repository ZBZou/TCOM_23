%% Ber DPC
clear all;
close all;
%% load
c = load ('h_v2x.mat');
par = load('par_v2x.mat');
%% para
Nu = par.Nu;
Nt = par.Nt;
Bw = par.Bw;
N_carrier = par.Ns;
fs = par.fs;
dt = par.dt;
T_diff = round(dt*Bw);
N_ofdm = 10;
N_cp = N_carrier/4;
Nu = size(c.h,1);
T = (N_carrier+N_cp)*N_ofdm;
SNRdB = [1:30];
N_loop = 5;
M = 16;
K = log2(M);

%% interpolate
ch = zeros(Nu,Nt,N_carrier,T);
for u1 = 1:Nu
   for u2 = 1:Nt
      for t1 = 1:N_carrier
         for t2 = 1:T
            t_temp = floor(t2/T_diff)+1;
            ch(u1,u2,t1,t2) = c.h(u1,u2,t1,t_temp);
         end
      end
   end
end
%%
Hf = fft(ch, [], 3);

%% CSI 
% H = c.h(1:Nu,1:Nu,:,1:T);

Hc = permute(ch, [1,4,2,3]);

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
%% sync
% for i = 1:Nu
%     [~,id] = max(Hst(i,1,i,:));
%     for j = 1:Nu
%         for t1 = 1:T
%             for t2 = t1:T-id+1
%                 Hst(i,t1,j,t2) = Hst(i,t1,j,t2+id-1);
%             end
%         end
%     end
% end
%% data
% data = randi([0 1],Nu,N_carrier*N_ofdm*K);
% s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
% s = s';
% s_temp = reshape(s, Nu, N_carrier, []);
% s_temp = permute(s_temp, [2,3,1]); %NFFT by T by Nu
% %% OFDM
% s_ofdm = ofdmmod(s_temp, N_carrier, N_cp);
% s_ofdm = permute(s_ofdm, [2,1]);
% %% DPC
% L = zeros(Nu,Nu,T);
% Q = zeros(Nu,Nu,T);
% Tx_signal = s_ofdm;
% xp = s_ofdm;
% %% 
% for i = 1:T
%    H = squeeze(Hst(:,i,:,i));
%    [Q_temp,R_temp] = qr(H');
%    L(:,:,i)=R_temp'; Q(:,:,i)=Q_temp';
%    % pre_eq
% %    pre_s = s_ofdm;
% %    pre_s(:,i) = inv(diag(diag(L(:,:,i))))*squeeze(s_ofdm(:,i));  
%    for m=2:Nu % Eqs.(13.39)(13.41)
%        xp(m,i) = xp(m,i) - L(m,1:m-1,i)/L(m,m,i)*squeeze(xp(1:m-1,i));
%    end
%    Tx_signal(:,i) = Q(:,:,i)'*xp(:,i);
%    
% end
% 
% %% Receiver
% Rx_signal = zeros(Nu,T);
% r = zeros(Nu,T);
% 
% % 
% for i = 1:Nu
%    for j = 1:T
%       temp_H = conj(squeeze(Hst(i,:,:,j))');
%       Rx_signal(i,j) = sum(sum(temp_H .* Tx_signal)); %+ sigma*(randn()+randn()*1i)/sqrt(2);
%    end
% end
% 
% %% equa
% for i = 1:T
%    in_L = inv(diag(diag(L(:,:,i))));
%    r(:,i) = in_L*Rx_signal(:,i);
% end
% %% OFDMdemod
% r_temp = permute(r, [2,1]);
% s_demod = ofdmdemod(r_temp, N_carrier, N_cp); % NFFT by T by Nu
% s_hat = reshape(permute(s_demod, [3,1,2]),Nu,[]);
% %% BER
% Ber_dpc = myber(data,s_hat,M)
BER = zeros(1,length(SNRdB));
for SNR = 1:length(SNRdB)
   sigma2 = 0.5*10^(-SNRdB(SNR)/10); sigma = sqrt(sigma2);
   for loop = 1:N_loop
      %% data
      data = randi([0 1],Nu,N_carrier*N_ofdm*K);
      s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
      s = s';
      %%
      s = reshape(s, Nu, N_carrier, []);
      %% DPC
         L = zeros(Nu,Nu,N_carrier, N_ofdm);
         Q = zeros(Nu,Nu,N_carrier, N_ofdm);
         Tx_signal = s;
         xp = s;
      for i = 1:N_ofdm 
         t = (i-1)*(N_carrier+N_cp)+1;
         for j = 1:N_carrier
            H = squeeze(Hf(:,:,j,t));
            [Q_temp,R_temp] = qr(H');
            L(:,:,j,i)=R_temp'; Q(:,:,j,i)=Q_temp';
            for m=2:Nu % Eqs.(13.39)(13.41)
               xp(m,j,i) = xp(m,j,i) - L(m,1:m-1,j,i)/L(m,m,j,i)*squeeze(xp(1:m-1,j,i));
            end
            Tx_signal(:,j,i) = Q(:,:,j,i)'*squeeze(xp(:,j,i));
         end
      end

      %% OFDM
      s_temp = permute(Tx_signal, [2,3,1]); %NFFT by T by Nu
      s_ofdm = ofdmmod(s_temp, N_carrier, N_cp);
      s_ofdm = permute(s_ofdm, [2,1]);

      %% Receiver
      % Rx_signal = zeros(Nu,N_carrier,N_ofdm);
      % r = zeros(Nu,N_carrier,N_ofdm);
      % for i = 1:N_carrier
      %    for j = 1:N_ofdm
      %       temp_t = (j-1)*(N_carrier)+1;
      %       temp_hf = squeeze(Hf(:,:,i,t));
      %       Rx_signal(:,i,j) = temp_hf*Tx_signal(:,i,j);
      %    end
      % end
      Rx_signal = zeros(Nu,T);
      r = zeros(Nu,N_carrier,N_ofdm);
      for i = 1:Nu
         for j = 1:T
            temp_H = conj(squeeze(Hst(i,:,:,j))');
            Rx_signal(i,j) = sum(sum(temp_H .* s_ofdm));
         end
      end

      %% OFDMdemod
      r_temp = permute(Rx_signal, [2,1]);
      s_demod = ofdmdemod(r_temp, N_carrier, N_cp); % NFFT by T by Nu
      % temp_s = reshape(permute(s_demod, [3,1,2]),Nu,[]);
      temp_s = permute(s_demod, [3,1,2]);
      %% demod
      % temp_s = Rx_signal;
      for u = 1:Nu
         for i = 1:N_ofdm
            for j = 1:N_carrier  
               temp_t = (i-1)*(N_carrier)+1;
               in_L = inv(diag(diag(L(u,u,j,i))));
               r(u,j,i) = in_L*temp_s(u,j,i);
            end
         end
      end

      s_hat = reshape(r, Nu, []);
      %% BER
      BER(SNR)=BER(SNR)+myber(data,s_hat,M);
   end
   BER(SNR) = BER(SNR)/N_loop;
end

savefile = 0;
if savefile == 1
   save DPC_ber.mat BER SNRdB;
end