%% plot cdfs
clear all;
close all;
%% load
c = load ('h_d50_v120.mat');
par = load('par_d50_v120.mat');
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
N_cp = N_carrier/4;
Nu = size(c.h,1);
T = (N_carrier+N_cp)*N_ofdm;
SNRdB = 30;
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
s_ofdm_hogmt = [s_ofdm_hogmt,s_ofdm_hogmt(:,1:1/4*size(s_ofdm_hogmt,2))];
%%
% [ecdf,pcdf,scdf] = eigencdfs(s_ofdm_hogmt, Kh);

s = s_ofdm_hogmt;
h = Kh;
Nu = size(h,1);
T = size(h,2);
H = reshape(h, Nu*T, Nu*T);
S = reshape(s,Nu*T,1);
X = zeros(Nu*T,1);
display("Decomposing Channel");
[U,temp_sig,V] = svd(H);
display("Decomposed");
%%
sig = diag(temp_sig);
Nsig = length(sig);
ecdf = zeros(Nsig,1);
pcdf = zeros(Nsig,1);
scdf = zeros(Nsig,1);
eigen = sig.^2;
s_coeff = 10*log10(abs(conj(S')*U).^2);
p = 10*log10(abs(s_coeff./conj(sig')).^2);
for i = 1:Nsig
   ecdf(i,1) = sum(eigen(1:i));
   scdf(i,1) = sum(s_coeff(1:i));
   pcdf(i,1) = sum(p(1:i));
end
ecdf = ecdf/sum(eigen);
pcdf = pcdf/sum(p);
scdf = scdf/sum(s_coeff);
%%
x = [1:Nsig]/8000;
f = figure;
% yyaxis left
plot(x,ecdf,x,scdf,x,pcdf,'LineWidth',2);
xlabel('Normalized Number of Eigenfunctions');
ylabel('Normalized Cumulative Function');
% yyaxis right
% plot(x,sig/sum(sig),x,s_coeff/sum(s_coeff),x,p/sum(p));
legend('Eigenvalue \lambda','Cancelled Energy e^\prime','Consumed Energy e','Location','southeast');
set(gca,'fontsize',20);
grid on;
[h, wd, ht] = tightfig();
filename = 'cumulative';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f, name1);
exportgraphics(f, name2);