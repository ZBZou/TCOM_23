%%%%%%Non-WSSUS channel%%%%

clear all;
close all,
disp('Probed channel matrix loading...')
load 'Z_True_Complex.mat';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel's TF transfer function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LH = Z_True_Complex;
t=[0:size(LH,1)-1];
f=[0:size(LH,2)-1];

figure
mesh(f,t,(abs(LH)));
view(-5,35);
xlabel('f');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 0 0.25])
ylabel('t');
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [0 0 0.25])
% title('TF transfer function','fontsize',14,'fontweight','b');
grid off;
Ax = gca;Ax.ZAxis.Visible = 'off';
Ax.ZGrid = 'off';Ax.Color = 'none';
set(gca,'fontsize',30)
xtics=[0,63];
ytics=[0,1000];
set(gca,'XTick',xtics,'YTick',ytics)
[h, wd, ht] = tightfig();
print -opengl -dpdf -r600 ht_tau.pdf
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel's time-varying response and envelope %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = f;
H_tau = ifft(LH,[],2); 
figure
% subplot(2,1,1);
mesh(tau,t,abs(H_tau));
xlabel('\tau');
xlabh = get(gca,'XLabel');

ylabel('t');
ylabh = get(gca,'YLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 0 0.25])
set(ylabh,'Position',get(ylabh,'Position') + [0 0 0.25])
grid off;
Ax = gca;Ax.ZAxis.Visible = 'off';
Ax.ZGrid = 'off';Ax.Color = 'none';
xtics=[0,63];
ytics=[0,1000];
view(-5,35);
set(gca,'XTick',xtics,'YTick',ytics)
% title('Time-varying response','fontsize',14,'fontweight','b');
set(gca,'fontsize',30)
% subplot(2,1,2);
% plot(t,abs(sum(H_tau,2)));
% grid on;
% xlabel('t (s)','fontsize',12);
% axis([0 60 0 3]);
% title('Envelope','fontsize',14,'fontweight','b');
[h, wd, ht] = tightfig();
print -opengl -dpdf -r600 ht_tau.pdf

%%
%plot the channel's spreading function %
nu=t;
SH = fft(H_tau,[],1);
figure
mesh((abs(SH)));
xlabel('\tau');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('\nu');
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
grid off;
Ax = gca;Ax.ZAxis.Visible = 'off';
Ax.ZGrid = 'off';Ax.Color = 'none';
set(gca,'fontsize',30)
xtics=[0,63];
ytics=[0,1000];
set(gca,'XTick',xtics,'YTick',ytics)
% title('Spreading function','fontsize',14,'fontweight','b');
%%
%Kernel
T1 = 100;
F1 = 64;
T2 = 100;
F2 = 64;
kH = Kernel(LH, T1, F1, T2, F2);
test_k = squeeze(kH(10,64,:,:));
figure
mesh(10*log10(abs(test_k)));


