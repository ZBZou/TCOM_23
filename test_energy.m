%%Energy compare
clear all;
close all;
clc;

Nu = 10;
M = 16;
K = log2(M);
data = randi([0 1],Nu,K);
s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
s = s';
ep = 0;
P_H = 1;
P_tx = Nu;
Nt = 1;
sq2 = sqrt(2);
SNRdB = [0:1:30];
N_loop = 1000;

BER_DPC = zeros(length(SNRdB),1);
BER_HOGMT = zeros(length(SNRdB),1);
P_DPC = zeros(length(SNRdB),1);
P_HOGMT = zeros(length(SNRdB),1);
%%
for snr = 1:length(SNRdB)
   sigma2 = 0.5*Nt*10^(-SNRdB(snr)/10); sigma = sqrt(sigma2);
   for n = 1:N_loop
      %% H generate
      H = zeros(Nu,Nu);
      for i = 1:Nu
         for j = 1:Nu
            H(i,j) = sqrt(P_H)*(randn()+randn()*1i)/sq2;
         end
      end
      %% Precoding
      X_HOGMT = HOGMT(s, H, ep);
      [X_DPC,L] = DPC(s,H);
      
      %% power constriant
      p_HOGMT = sum(abs(X_HOGMT).^2);
      p_DPC = sum(abs(X_DPC).^2);
      X_HOGMT = X_HOGMT * sqrt(P_tx/p_HOGMT);
      X_DPC = X_DPC * sqrt(P_tx/p_DPC);
      %% Rx
      r_HOGMT = H*X_HOGMT + sigma*(randn()+randn()*1i)/sq2;
      r_DPC = H*X_DPC + sigma*(randn()+randn()*1i)/sq2;

      %% Normalization at Rx
      r_HOGMT = r_HOGMT * sqrt(p_HOGMT/P_tx);
      r_DPC = r_DPC * sqrt(p_DPC/P_tx);
      %% decoding for DPC
      r_DPC = inv(diag(diag(L)))*r_DPC;
      
      %% BER
      ne_HOGMT = myber(data,r_HOGMT,M);
      ne_DPC = myber(data,r_DPC,M);

      ber_HOGMT = ne_HOGMT/(Nu*K);
      ber_DPC = ne_DPC/(Nu*K);
      
      %% power constriant
      p_HOGMT = sum(abs(X_HOGMT).^2);
      p_DPC = sum(abs(X_DPC).^2);
      
      %% loop result 
      BER_DPC(snr) = BER_DPC(snr)+ber_DPC;
      BER_HOGMT(snr) = BER_HOGMT(snr)+ber_HOGMT;
      P_DPC(snr) = P_DPC(snr)+p_DPC;
      P_HOGMT(snr) = P_HOGMT(snr)+p_HOGMT;
   end
   BER_DPC(snr) = BER_DPC(snr)/N_loop;
   BER_HOGMT(snr) = BER_HOGMT(snr)/N_loop;
   P_DPC(snr) = P_DPC(snr)/(N_loop*Nu);
   P_HOGMT(snr) = P_HOGMT(snr)/(N_loop*Nu);
end

%% plt BER
figure;
semilogy(SNRdB, BER_DPC, SNRdB, BER_HOGMT);
legend('DPC','HOGMT-Precoding');
title('BER');
%% plt power
figure;
semilogy(SNRdB, P_DPC, SNRdB, P_HOGMT);
title('Power')
legend('DPC','HOGMT-Precoding');

%% plt channel
figure;
surf(abs(H));
title('Channel H');

%% one time HOGMT 

SNRdB = 20;
sigma2 = 0.5*Nt*10^(-SNRdB/10);
sigma = sqrt(sigma2);

H = zeros(Nu,Nu);
for i = 1:Nu
   for j = 1:Nu
      H(i,j) = sqrt(P_H)*(randn()+randn()*1i)/sq2;
   end
end

h_r = real(H);
h_i = imag(H);
h = [h_r, -h_i;h_i, h_r];
[psi,temp_sig,phi] = svd(h);

Sig = diag(temp_sig);
N = length(Sig); 
sn = [];
xn = [];

sr = real(s);
si = imag(s);
S = [sr;si];

X = zeros(length(S),1);

for i = 1:N
   sn(i) = dot(S, psi(:,i));
   xn(i) = sn(i)/Sig(i);
   X = X + xn(i)*phi(:,i);
end

Psi_r = psi(1:N/2,:);
Psi_i = psi(N/2+1:end,:)*1i;
Psi = Psi_r + Psi_i;

Phi_r = phi(1:N/2,:);
Phi_i = phi(N/2+1:end,:)*1i;
Phi = Phi_r + Phi_i;

Xr = X(1:N/2);
Xi = X(N/2+1:end)*1i;
X_HOGMT = Xr+Xi;

r_HOGMT = H*X_HOGMT + sigma*(randn()+randn()*1i)/sq2;
r_r = real(r_HOGMT);
r_i = imag(r_HOGMT);

%% plot coefficients
figure;
subplot(3,1,1);
plot([1:N], Sig.^2);
title('\sigma_n^2');
xlim([1 N]);
subplot(3,1,2);
plot([1:N], xn.^2);
title('x_n^2');
xlim([1 N]);
subplot(3,1,3);
plot([1:N], sn.^2);
title('s_n^2');
xlim([1 N]);

%% plot eigenfunctions
figure;
subplot(2,1,1)
[x1,y1] = meshgrid(1:size(Phi,2), 1:size(Phi,1));
plot3(x1,y1, abs(Phi));
view(-8,67);
xlabel('n');
xlim([1 N]);
title('|\Phi|');
subplot(2,1,2)
[x2,y2] = meshgrid(1:size(Psi,2), 1:size(Psi,1));
plot3(x2,y2, abs(Psi));
view(-8,67);
xlim([1 N]);
xlabel('n');
title('|\Psi|');

%% plot decomposition of X and S
figure;
subplot(2,1,1)
X_decomp = [X_HOGMT,Phi.*xn];
[x1,y1] = meshgrid(0:size(X_decomp,2)-1, 1:size(X_decomp,1));
plot3(x1,y1, abs(X_decomp));
view(-8,67);
ylabel('u')
xlabel('n');
xlim([0 N]);
str1 ='$$ |X = \sum_n x_n \phi_n| $$';
title(str1,'Interpreter','latex')
subplot(2,1,2)
s_decomp = [s,Psi.*sn];
[x2,y2] = meshgrid(0:size(s_decomp,2)-1, 1:size(s_decomp,1));
plot3(x2,y2, abs(s_decomp));
view(-8,67);
xlim([0 N]);
ylabel('u')
xlabel('n');
str2 ='$$ |S = \sum_n s_n \psi_n| $$';
title(str2,'Interpreter','latex')

%% optimize ber (intererence ratio and power ratio)
Pow = sum(xn.^2);
pr = [];
pi = [];
for i = 1:N
   pr(i) = xn(i).^2/(Pow-xn(i));
end

var_int = sn.^2/Nu;
SNRdB_temp = [5];
var_improve = zeros(N, length(SNRdB_temp));
var_diff = zeros(N, length(SNRdB_temp));

for i = 1:length(SNRdB_temp)
   var = 0.5*Nt*10^(-SNRdB_temp(i)/10);
   for j = 1:N
       var_noise(j,i) = var*pr(j);
   end
   var_diff(:,i) = var_noise(:,i) - var_int';
end

figure;
plot([1:N], var_int, [1:N], var_noise(:,1));
legend('interference','noise')
xlabel('n')
ylabel('var')