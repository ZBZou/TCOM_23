%% statistics 
close all
clear all
%%
c = load ('h_d10_v120.mat');
par = load('par_d10_v120.mat');

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
hc = c.h;
T_diff = round(dt*Bw);
N_ofdm = 10;
N_cp = Ns/4;
T = (Ns+N_cp)*N_ofdm; % 10 ofdm symbols = 0.04 ms
% Stat_T = 100*T;
%% interpolate
% hc = zeros(Nu,Nt,N_carrier,T);
% for u1 = 1:Nu
%    for u2 = 1:Nt
%       for t1 = 1:N_carrier
%          for t2 = 1:Stat_T
%             t_temp = round(t2/T_diff)+1;
%             hc(u1,u2,t1,t2) = h(u1,u2,t1,t_temp);
%          end
%       end
%    end
% end
%%
t1 = round(1e-3/dt);
plt_h_t = abs(squeeze(hc(1,1,:,:)));
acf = [];
lags = [];
t = [1:t1];
for i = 1:length(t)
    [acf(:,i),lags(:,i)] = autocorr(plt_h_t(:,t(i)));
end
acf = permute(acf,[2,1]);
[xcor, ycor] = meshgrid(1:size(acf,1), 1:size(acf,2));
figure;
% plot3(ycor, xcor/2, acf,'-o');
mesh(ycor, xcor/t1, acf');
% legend('1ms','2ms','3ms','4ms','5ms','6ms','7ms','8ms','9ms','10ms');
% plot(lags, acf(1,:),lags, acf(2,:),lags, acf(3,:),lags, acf(4,:),lags, acf(5,:));
xlabel('Lags');
ylabel('Time (ms)');
zlabel('ACF');
view(35,35);
set(gca,'fontsize',16);
filename = 'Acf';
savefile = 0;
[h, wd, ht] = tightfig();
if savefile == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%% boxplot
figure;
boxplot(acf);
xlabel('Lags');
ylabel('Average ACF');
set(gca,'fontsize',16);
filename = 'boxplot';
savefile = 0;
[h, wd, ht] = tightfig();
if savefile == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%% CMD
% H_temp = sum(hc,3);
% plot_dtx = [];
% plot_drx = [];
% for j = 1:t1
% k1 = 1;
% k2 = j;
% L = t1;
% Rtx1 = zeros(10,10);
% Rrx1 = zeros(10,10); 
% Rtx2 = zeros(10,10); 
% Rrx2 = zeros(10,10); 
% for i = 1:L
%    Hk1 = squeeze(H_temp(:,:,k1));
%    Hk2 = squeeze(H_temp(:,:,k2));
%    Rtx1 = Rtx1 + conj(Hk1')*conj(Hk1);
%    Rrx1 = Rrx1 + Hk1*Hk1';
%    Rtx2 = Rtx2 + conj(Hk2')*conj(Hk2);
%    Rrx2 = Rrx2 + Hk2*Hk2';
% end
% Rtx1 = Rtx1/L; 
% Rrx1 = Rrx1/L; 
% Rtx2 = Rtx2/L;
% Rrx2 = Rrx2/L; 
% dtx = abs(1-(trace(Rtx1*Rtx2)/(norm(Rtx1,"fro")*norm(Rtx2,"fro")))); 
% drx = abs(1-(trace(Rrx1*Rrx2)/(norm(Rrx1,"fro")*norm(Rrx2,"fro")))); 
% plot_dtx = [plot_dtx;dtx];
% plot_drx = [plot_drx;drx];
% end
% %%
% xl = [1:t1]*dt*1e6;
% plot(xl,plot_dtx,xl,plot_drx);
% xlabel('Time (\mu s)');
% ylabel('Correlation Matrix Distance');
% legend('Tx','Rx','location','southeast');
% xlim([0,1000]);
% set(gca,'fontsize',20);
t1 = 5*round(1e-3/dt);

H_temp = sum(hc,3);
plot_dtx = [];
plot_drx = [];
for k1 = t1+1:2*t1
   for k2 = k1-t1:k1+t1-1
   L = t1;
   Rtx1 = zeros(10,10);
   Rrx1 = zeros(10,10); 
   Rtx2 = zeros(10,10); 
   Rrx2 = zeros(10,10); 
   for i = 1:L
      Hk1 = squeeze(H_temp(:,:,k1));
      Hk2 = squeeze(H_temp(:,:,k2));
      Rtx1 = Rtx1 + conj(Hk1')*conj(Hk1);
      Rrx1 = Rrx1 + Hk1*Hk1';
      Rtx2 = Rtx2 + conj(Hk2')*conj(Hk2);
      Rrx2 = Rrx2 + Hk2*Hk2';
   end
   Rtx1 = Rtx1/L; 
   Rrx1 = Rrx1/L; 
   Rtx2 = Rtx2/L;
   Rrx2 = Rrx2/L; 
   dtx = abs(1-(trace(Rtx1*Rtx2)/(norm(Rtx1,"fro")*norm(Rtx2,"fro")))); 
   drx = abs(1-(trace(Rrx1*Rrx2)/(norm(Rrx1,"fro")*norm(Rrx2,"fro")))); 
   plot_dtx(k1,k2) = dtx;
   plot_drx(k1,k2) = drx;
   end
end

%%
figure;
plot_temp = zeros(t1, 2*t1);
for i = 1:t1
   for j = 1:2*t1
      plot_temp(i,j) = plot_dtx(t1+i, i+j-1);
   end
end
imagesc(plot_temp');colorbar;
colormap('hot');
xtick = [1:size(plot_temp,1)/4:size(plot_temp,1)];
xlab = round([0 250 500 750 1000]);
ytick = [1:size(plot_temp,2)/8:size(plot_temp,2)];
ylab = round([-1000 -750 -500 -250 0 250 500 750 1000]);

set(gca,'XTick',xtick);
set(gca,'XTickLabel',xlab);
set(gca,'YTick',ytick);
set(gca,'YTickLabel',ylab);
xlabel('t (\mu s)');
ylabel('\Delta t (\mu s)');
zlabel('Correlation Matrix Distance')
set(gca,'fontsize',26);
hold on;
[C,h1] = contour(plot_temp',[0.2,0.7],'b--','ShowText','on','LineWidth',2);
clabel(C,h1,'FontSize',20,'Color','blue');
filename = 'CMD_TX';
savefile = 0;
if savefile == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end
%%
figure;
plot_temp = zeros(t1, 2*t1);
for i = 1:t1
   for j = 1:2*t1
      plot_temp(i,j) = plot_drx(t1+i, i+j-1);
   end
end
imagesc(plot_temp');colorbar;
colormap('hot');
xtick = [1:size(plot_temp,1)/4:size(plot_temp,1)];
xlab = round([0 250 500 750 1000]);
ytick = [1:size(plot_temp,2)/8:size(plot_temp,2)];
ylab = round([-1000 -750 -500 -250 0 250 500 750 1000]);

set(gca,'XTick',xtick);
set(gca,'XTickLabel',xlab);
set(gca,'YTick',ytick);
set(gca,'YTickLabel',ylab);
xlabel('t (\mu s)');
ylabel('\Delta t (\mu s)');
zlabel('Correlation Matrix Distance')
set(gca,'fontsize',26);
hold on;
[C,h1] = contour(plot_temp',[0.2,0.3],'b--','ShowText','on','LineWidth',2);
clabel(C,h1,'FontSize',20,'Color','blue');
filename = 'CMD_RX';
savefile = 0;
if savefile == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gca, name1);
    exportgraphics(gca, name2);
end