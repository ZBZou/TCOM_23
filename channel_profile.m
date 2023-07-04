%% channel profile
clear all;
close all;

%%
c = load ('h_d50_v120.mat');
par = load('par_d50_v120.mat');
%%
Nu = par.Nu;
Nt = par.Nt;
Bw = par.Bw;
Ns = par.Ns;
fs = par.fs;
dt = par.dt;
save = 1;
N_snap = 4000;
fc = par.fc;
p = par.p;
h = c.h;
N_carrier = par.Ns;
N_ofdm = 10;
N_cp = N_carrier/4;
T = (N_carrier+N_cp)*N_ofdm;
T_diff = round(dt*Bw);
%% clusters
p(1).visualize_clusters([],[],1);
set(gca,'fontsize',14);
filename = 'cluster';
if save == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    print(gcf, filename, '-dpdf', '-bestfit')
end
%% ht
save_ht = 0;
ht_name = 'ht_u1_ant1';
typet = 'time-delay';
plt_h_t = squeeze(h(1,1,:,1:300));
plot_h(plt_h_t, dt*1e3, 1e6/Bw, typet, save_ht, ht_name);
%% plt hs
save_hs = 0;
hs_name = 'hs_1s';
typet = 'space-space';
plt_h_s = squeeze(h(:,:,1,1));
plot_h(plt_h_s, 1, 1, typet, save_hs, hs_name);
%% plt hst
save_hst = 0;
hst_name = 'hst_2u_1ms';
typef = 'space-time';
t1 = round(1e-3/dt);
t2 = 2*round(1e-3/dt);
t3 = 100*round(1e-3/dt);
t4 = 1000*round(1e-3/dt);
plt_h_st = squeeze(h(2,:,:,t2));
plot_h(plt_h_st, 1e6/Bw, 1, typef, save_hst, hst_name);
%% plt hf
hf = fft(h,[],3);
save_ht = 0;
ht_name = 'example_hf_v2x';
typef = 'time-frequency';
plt_h_f = squeeze(hf(1,1,:,1:200));
plt_h_f = plt_h_f;
f_cor = (fc-Bw/2+[1:size(plt_h_f,1)]*Bw/size(plt_h_f,1))/1e6;
plot_h(plt_h_f, dt*1e3, f_cor, typef, save_ht, ht_name);
%% plt kernel
save_k = 1;
hk_name = 'hst_1u_1ms';
typef = 'kernel';
t1 = round(1e-3/dt);
t2 = 2*round(1e-3/dt);
t3 = 100*round(1e-3/dt);
t4 = 1000*round(1e-3/dt);
plt_k = squeeze(h(1,:,:,t1));
plot_h(plt_k, 1e6/Bw, 1, typef, save_k, hk_name);
%% plt pdp
save_pdp = 0;
plt_pdp = abs(squeeze(h(1,1,:,:))).^2;
plt_pdp = 10*log10(plt_pdp');
imagesc(plt_pdp);
colormap('hot');
colorbar;
set(gca,'YTick',1:15625:size(plt_pdp,1));
set(gca,'YTickLabel',(0:15625:size(plt_pdp,1))*dt*1e3);
set(gca,'XTick',1:10:size(plt_pdp,2));
set(gca,'XTickLabel',(0:10:size(plt_pdp,2))*1e6/Bw);
xlabel('Delay (\mus)');
ylabel('Time (ms)');
set(gca,'fontsize',18);
filename = 'pdp';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(gcf, name1);
exportgraphics(gcf, name2);
%% plt Doppler
% save_ds = 0;
% hst_name = 'ds_v2x';
% typef = 'ds';
% temp_ds = squeeze(h(1,1,:,:));
% temp_ds = ifft2(temp_ds);
% temp_ds = fftshift(temp_ds);
% plt_ds = abs(temp_ds).^2;
% plt_ds = 10*log10(plt_ds');
% imagesc(plt_ds);
% colormap('hot');
% colorbar;
% % set(gca,'YTick',1:15625:size(plt_pdp,1));
% % set(gca,'YTickLabel',(0:15625:size(plt_pdp,1))*dt*1e3);
% % set(gca,'XTick',1:10:size(plt_pdp,2));
% % set(gca,'XTickLabel',(0:10:size(plt_pdp,2))*1e6/Bw);
% xlabel('Delay (\mus)');
% ylabel('Time (ms)');
% set(gca,'fontsize',16);
%% plt kernel

hc = zeros(Nu,Nt,N_carrier,T);
for u1 = 1:Nu
   for u2 = 1:Nt
      for t1 = 1:N_carrier
         for t2 = 1:T
            t_temp = round(t2/T_diff)+1;
            hc(u1,u2,t1,t2) = h(u1,u2,t1,t_temp);
         end
      end
   end
end

% kernl 

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


save_kt = 0;
hst_name = 'kt_u1_ant1';
typeK = 'space-time';
t1 = round(1e-3/dt);
t2 = 2*round(1e-3/dt);
t3 = 100*round(1e-3/dt);
t4 = 1000*round(1e-3/dt);

plt_k_t = squeeze(Hst(1,:,:,t2));
plot_h(plt_k_t, dt*1e3, dt*1e3, typeK, save_kt, kt_name);

%% plt eigen functions
% interpolate
hc = zeros(Nu,Nt,N_carrier,T);
for u1 = 1:Nu
   for u2 = 1:Nt
      for t1 = 1:N_carrier
         for t2 = 1:T
            t_temp = round(t2/T_diff)+1;
            hc(u1,u2,t1,t2) = h(u1,u2,t1,t_temp);
         end
      end
   end
end

% kernl 

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
% 
% save_kt = 0;
% kt_name = 'kt_u1_ant1';
% typet = 'time-time';
% plt_k_t = squeeze(Hst(1,1:100,1,1:100));
% plot_h(plt_k_t, dt*1e3, dt*1e3, typet, save_kt, kt_name);

%% map

H = reshape(Hst, Nu*T, Nu*T);

%% svd

display("Decomposing Channel");
[temp_psi,temp_sig,temp_phi] = svd(H);
display("Decomposed");

%% mapback
phi_1 = abs(reshape(temp_psi(:,1),Nu,T));
phi_2 = abs(reshape(temp_psi(:,2),Nu,T));
psi_1 = abs(reshape(temp_phi(:,1),Nu,T));
psi_2 = abs(reshape(temp_phi(:,2),Nu,T));
%% plt
f1 = figure;
mesh(phi_1');
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(45,50);
set(gca,'fontsize',20)
% title('\phi_1 (u^\prime,t^\prime)');
% axis off;
filename = 'phi_1';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f1, name1);
exportgraphics(f1, name2);
%%
f2 = figure;
mesh(phi_2');
grid on;
xlabel('u^\prime');
ylabel('t^\prime');
view(45,50);
set(gca,'fontsize',20);
% title('\phi_2 (u^\prime,t^\prime)');
% axis off;
filename = 'phi_2';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f2, name1);
exportgraphics(f2, name2);
%%
f3 = figure;
mesh(psi_1');
grid on;
xlabel('u');
ylabel('t');
view(45,50);
set(gca,'fontsize',20)
% title('\psi_1 (u,t)');
% axis off;
filename = 'psi_1';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f3, name1);
exportgraphics(f3, name2);
%%
f4 = figure;
mesh(psi_2');
grid on;
xlabel('u');
ylabel('t');
view(45,50);
set(gca,'fontsize',20);
% title('\psi_2 (u,t)');
% axis off;

filename = 'psi_2';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f4, name1);
exportgraphics(f4, name2);
%%
% figure;
% subplot(2,2,1)
% mesh(phi_1);
% grid on;
% xlabel('u^\prime');
% ylabel('t^\prime');
% view(-15,50);
% set(gca,'fontsize',20)
% title('\phi_1 (u^\prime,t^\prime)');
% % axis off;
% subplot(2,2,3)
% mesh(phi_2);
% grid on;
% xlabel('u^\prime');
% ylabel('t^\prime');
% view(-15,50);
% title('\phi_2 (u^\prime,t^\prime)');
% set(gca,'fontsize',20)
% % axis off;
% subplot(2,2,2)
% mesh(psi_1);
% grid on;
% xlabel('u');
% ylabel('t');
% view(-15,50);
% set(gca,'fontsize',20)
% title('\psi_1 (u,t)');
% % axis off;
% subplot(2,2,4)
% mesh(psi_2);
% grid on;
% xlabel('u');
% ylabel('t');
% title('\psi_2 (u,t)');
% view(-15,50);
% set(gca,'fontsize',20)
% % axis off;

%% power map

