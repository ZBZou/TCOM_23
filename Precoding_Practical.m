%% Precoding for Practical Scenarios
close all
clear all

%% load
c = load ('h_v2x.mat');
par = load('par_v2x.mat');
%% parameters
Nu = par.Nu;
Nt = par.Nt;
Bw = par.Bw;
N_carrier = par.Ns;
fs = par.fs;
dt = par.dt;
T_diff = round(dt*Bw);
% Bw = 20e6;
% N_carrier = 64;
N_ofdm = 10;
r_cp = 0;
N_cp = N_carrier*r_cp;
Nu = size(c.h,1);
T = (N_carrier+N_cp)*N_ofdm;
SNRdB = 25;
M = 16;
K = log2(M);
sigPow = 0;
sigma2 = 0.5*10^(-SNRdB/10); sigma = sqrt(sigma2);
%% interpolate
H = zeros(Nu,Nt,N_carrier,T);
for u1 = 1:Nu
   for u2 = 1:Nt
      for t1 = 1:N_carrier
         for t2 = 1:T
            t_temp = round(t2/T_diff)+1;
            H(u1,u2,t1,t2) = c.h(u1,u2,t1,t_temp);
         end
      end
   end
end

%% CSI 
% H = c.h(1:Nu,1:Nu,:,1:T);

H = permute(H, [1,4,2,3]);

Kh = zeros(Nu,T,Nu,T);

for u1 = 1:Nu
    for u2 = 1:Nu
        for t1 = 1:T
            for t2 = t1:t1+size(H,4)-1
                if t2 > T
                    continue;
                end
                Kh(u1,t1,u2,t2) = H(u1,t1,u2,t2-t1+1);
            end
        end
    end
end
%% data
data = randi([0 1],Nu,N_carrier*N_ofdm*K);
s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
s = s';
%% OFDM
s_temp = reshape(s, Nu, N_carrier, []);
s_temp = permute(s_temp, [2,3,1]);
s_ofdm_temp = ofdmmod(s_temp, N_carrier, N_cp);
s_ofdm_temp = permute(s_ofdm_temp, [2,1]);
%%
s_ofdm_hogmt = ofdmmod(s_temp, N_carrier, 0);
s_ofdm_hogmt = permute(s_ofdm_hogmt, [2,1]);
%% 
s_ofdm_hogmt = [s_ofdm_hogmt,s_ofdm_hogmt(:,1:r_cp*size(s_ofdm_hogmt,2))];
%% HOGMT
X_hogmt = HOGMT(s_ofdm_hogmt, Kh, 1);

%% receiver
Rx = zeros(Nu, T);
for i = 1:Nu
   for j = 1:T
      Rx(i,j) = sum(sum(squeeze(Kh(i,j,:,:)) .* X_hogmt));% + sigma*(randn()+randn()*1i)/sqrt(2);
   end
end 

%% AWGN
sigPow = sum(sum(abs(Rx)))/(size(Rx,1)*size(Rx,2));
for i = 1:Nu
   for j = 1:T
      Rx(i,j) = Rx(i,j) + sigPow*sigma*(randn()+randn()*1i)/sqrt(2);
   end
end 
%% demod
temp_s = Rx(:,1:1/(1+r_cp)*size(Rx,2));
r_temp = permute(temp_s, [2,1]);
s_demod = ofdmdemod(r_temp, N_carrier, 0);
s_hat = reshape(permute(s_demod, [3,1,2]),Nu,[]);
ber = myber(data,s_hat,M)
% IdeaBer(data, s_ofdm_temp, M, Nu, T, sigma, N_carrier, N_cp)