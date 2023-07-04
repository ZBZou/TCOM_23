%%%%%%Non-WSSUS channel%%%%

clear all;
close all,
disp('Probed channel matrix loading...')
load 'channel_demo.mat';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel's time-varying response and envelope %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau=[0:size(g_lk,1)-1]*Ts2;
t=[0:size(g_lk,2)-1]*Ts1;
figure
subplot(2,1,1);
mesh(tau*1E3,t,abs(g_lk).');
view(60,20);
axis([0 6 0 60 0 1]);
xlabel('\tau (ms)','fontsize',12);
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t (s)','fontsize',12);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
title('Time-varying response','fontsize',14,'fontweight','b');
subplot(2,1,2);
plot(t,abs(sum(g_lk)));
grid on;
xlabel('t (s)','fontsize',12);
axis([0 60 0 3]);
title('Envelope','fontsize',14,'fontweight','b');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel's TF transfer function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=[0:size(g_lk,1)-1]/(Ts2*10e3) ;
t=[0:size(g_lk,2)-1]*Ts1;
LH = fft(g_lk,[],1);
figure
mesh(f*1E3,t,10*log10(abs(LH)).');
xlabel('f (MHz)','fontsize',12);
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t (s)','fontsize',12);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
title('TF transfer function','fontsize',14,'fontweight','b');
%%
%plot the channel's spreading function %
tau=[0:size(g_lk,1)-1]*Ts2;
nu=([0:size(g_lk,2)-1]-size(g_lk,2)/2)*Ts1;
SH = fft(g_lk,[],2);
figure
mesh(tau*1E3,nu,10*log10(abs(SH)).');
view(60,20);
xlabel('\tau (ms)','fontsize',12);
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('\nu (Hz)','fontsize',12);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
title('Spreading function','fontsize',14,'fontweight','b');
%%
%Atomic coefficient
T1 = 100;
F1 = 10;
T2 = 100;
F2 = 10;
H_G = atomic_coefficient(LH', T1, F1, T2, F2);
test_H = squeeze(H_G(1,1,:,:));
figure
mesh(10*log10(abs(test_H)));

