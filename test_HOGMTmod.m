% %% HOGMT VS SVD
% clear all;
% close all;
% 
% Nu = 10;
% M = 16;
% K = log2(M);
% data = randi([0 1],Nu*K,1);
% SNR = [0:1:20];
% N_loop = 1000;
% 
% P_tx = Nu;
% 
% Ber_svd = [];
% Ber_HOGMT = [];
% for snr = 1:length(SNR)
%    ne_svd = 0;
%    ne_HOGMT = 0;
%    sigma2 = 0.5*10^(-SNR(snr)/10); sigma = sqrt(sigma2);
%    for n = 1:N_loop
%       H = zeros(Nu,Nu);
%       for i = 1:Nu
%          for j = 1:Nu
%             H(i,j) = (randn()+randn()*1i)/sqrt(2);
%          end
%       end
% 
%       %% SVD
% 
%       s_qam = qammod(reshape(data,[],Nu),M,'InputType','bit','UnitAveragePower',true);
%       s_qam = s_qam';
% 
%       [U,temp_sig1,V] = svd(H);
%       Sig1 = diag(temp_sig1);
%       X_svd = V*s_qam;
% 
%       %% HOGMT
%       data_temp = reshape(data,2*Nu,[]);
%       L = size(data_temp,2);
%       s_temp = zeros(2*Nu,1);
%       for i = 1:2*Nu
%          for j = 1:L
%             s_temp(i) = s_temp(i) + 2^(L-j)*data_temp(i,j);
%          end
%       end
% 
%       s_pam = pammod(s_temp,sqrt(M),0);
% 
% 
%       h_r = real(H);
%       h_i = imag(H);
%       h = [h_r, -h_i;h_i, h_r];
%       [psi,temp_sig2,phi] = svd(h);
%       Sig2 = diag(temp_sig2);
% 
%       N = size(psi,1);
%       Psi_r = psi(1:N/2,:);
%       Psi_i = psi(N/2+1:end,:)*1i;
%       Psi = Psi_r + Psi_i;
% 
% 
%       Phi_r = phi(1:N/2,:);
%       Phi_i = phi(N/2+1:end,:)*1i;
%       Phi = Phi_r + Phi_i;
% 
%       temp_s = reshape(s_pam,[],1);
% 
%       x = Psi*temp_s;
%       xr = real(x);
%       xi = imag(x);
% 
%       X = [xr;xi];
% 
% 
%       X_HOGMT = Phi*temp_s;
%       %% power constraint
%       p_HOGMT = sum(abs(X_HOGMT).^2);
%       p_SVD = sum(abs(X_svd).^2);
%       X_HOGMT = X_HOGMT * sqrt(P_tx/p_HOGMT);
%       X_svd = X_svd * sqrt(P_tx/p_SVD);
%       
%       %% rx
%       r_svd = H*X_svd + sigma*(randn()+1i*randn())/sqrt(2);
%       r_HOGMT = H*X_HOGMT + sigma*(randn()+1i*randn())/sqrt(2);
%       
%       %% Normalization at Rx
%       r_HOGMT = r_HOGMT * sqrt(p_HOGMT/P_tx);
%       r_svd = r_svd * sqrt(p_SVD/P_tx);
%       %% decode
%       y_svd = (U'*r_svd)./Sig1;
% 
%       r_r = real(r_HOGMT);
%       r_i = imag(r_HOGMT);
%       R = [r_r;r_i];
% 
%       y_HOGMT = psi'*R./Sig2;
% 
%       hat_s_svd = qamdemod(y_svd, M,'OutputType','bit','UnitAveragePower',true);
%       s_svd = qamdemod(s_qam, M,'OutputType','bit','UnitAveragePower',true);
% 
%       hat_s_HOGMT = pamdemod(y_HOGMT, sqrt(M), 0);
%       s_HOGMT = pamdemod(s_pam, sqrt(M), 0);
% 
%       ne_svd = ne_svd + length(find(reshape(s_svd,[],1) - reshape(hat_s_svd,[],1)));
%       ne_HOGMT = ne_HOGMT + length(find(reshape(s_HOGMT,[],1) - reshape(hat_s_HOGMT,[],1)));
%    end
%    Ber_svd(snr) = ne_svd/(N_loop*Nu*K);
%    Ber_HOGMT(snr) = ne_HOGMT/(N_loop*Nu*K);
% end
% 
% figure;
% semilogy(SNR, Ber_svd, SNR, Ber_HOGMT);
% legend('SVD','HOGMT')
% xlabel('snr')
% ylabel('ber')

%% H2

Nu = 2;
M = 4;
K = log2(M);
s = [1, 3; 2, 4];
% s = [1+1i, 3+3i; 2+2i, 4+4i];

H2 = zeros(Nu,Nu);
for i = 1:Nu
   for j = 1:Nu
      H2(i,j) = (randn()+randn()*1i)/sqrt(2);
   end
end

%% SVD

[U,temp_sig1,V] = svd(H2);
Sig1 = diag(temp_sig1);
X_svd = V*s;

%% HOGMT

h_r = real(H2);
h_i = imag(H2);
h = [h_r, -h_i;h_i, h_r];
[psi,temp_sig2,phi] = svd(h);
Sig2 = diag(temp_sig2);

N = length(Sig2);
Psi_r = psi(1:N/2,:);
Psi_i = psi(N/2+1:end,:)*1i;
Psi = Psi_r + Psi_i;


Phi_r = phi(1:N/2,:);
Phi_i = phi(N/2+1:end,:)*1i;
Phi = Phi_r + Phi_i;


temp_s = reshape(s,[],1);
X_HOGMT = Phi*temp_s;
%% rx
r_svd = H2*X_svd;
r_HOGMT = H2*X_HOGMT;

%% decode
y_svd = U'*r_svd./Sig1;

r_r = real(r_HOGMT);
r_i = imag(r_HOGMT);
R = [r_r;r_i];

y_HOGMT = Psi'*r_HOGMT./Sig2

