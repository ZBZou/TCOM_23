%%test time-varying response NS case
clear all;
close all,

T = 100; %time length
Taps = randi([10,20],1,T); %time-varying delay taps within [10,20]
T2 = 80; %data length
data = -3 + 2*randi([0,3],T2,1); %data serie
data = [data(T2-(T-T2)+1:end);data]; % data with cp

%NS channel generate
Stat = zeros(T,T);
H = zeros(T,T);
for i = 1:T
   vi = sqrt(i/T)*randn();
   ei = sqrt(i/T)*randn();
   H(i,1:Taps(i)) = ei + vi*randn(1,Taps(i));
   Stat(i,1:T) = ei + vi*randn(1,T);
end 

boxplot(Stat');
hold on
xlabel('Time index');
ylabel('Averaging channel Gains');
xtics=[1:10:T];
xticlab={'0' '10' '20' '30' '40' '50' '60' '70' '80' '90' '100'};
set(gca,'XTick',xtics,'XTickLabel',xticlab,'fontsize',30)



figure
mesh(10*log10(abs(H)));
xlabel('\tau','fontsize',14);
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t','fontsize',14);
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
set(gca,'fontsize',30)
% Ax = gca;Ax.ZAxis.Visible = 'off';
% Ax.ZGrid = 'off';Ax.Color = 'none';
% grid off;
% title('Time-varying response','fontsize',14,'fontweight','b');

% figure
% subplot(2,1,1);
% mesh(H);
% xlabel('\tau','fontsize',12);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
% ylabel('t','fontsize',12);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
% title('Time-varying response','fontsize',14,'fontweight','b');
% subplot(2,1,2);
% plot(10*log10(sum(H,2)));
% grid on;
% xlabel('t (s)','fontsize',12);
% title('Envelope','fontsize',14,'fontweight','b');

%% 2D kernel K(t,t'), t' = t-tau;
K = zeros(T,T);
for i = 1:T
   for j = 1:T
      if i-j < 0 || j > Taps(i)
         continue
      end
      K(i,i+1-j) = H(i,j);
   end
end

figure;
mesh(10*log10(abs(K)'));
xlabel('t^\prime');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t');
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
% title('2D Kernel K(t1,t2)','fontsize',14,'fontweight','b');
set(gca,'fontsize',30);

%% 2D Generalized Mercer's theorem
[U, S, V] = svd(K);
A = diag(S);
sig = [];
for i = 1:length(A)
   if A(i) < 0.01 
      N = i-1; % N most contributed eigenfunctions
      break;
   end
   sig(i) = A(i);
end 

phi = zeros(N,T);
psi = zeros(N,T);

for i = 1:N
   psi(i,:) = U(:,i);
   phi(i,:) = V(:,i);
end

figure;
subplot(2,1,1);
[x1,y1] = meshgrid(1:size(psi,2), 1:size(psi,1));  % T is your table
plot3(x1,y1,psi);
view(-8,67);
xlabel('N');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t');
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
set(gca,'fontsize',20);
title('Eigenfunctions \psi (t)','fontsize',14,'fontweight','b');
subplot(2,1,2);
[x2,y2] = meshgrid(1:size(phi,2), 1:size(phi,1));  % T is your table
plot3(x2,y2,phi);
view(-8,67);
xlabel('N');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('t^\prime');
ylabh = get(gca,'YLabel');
set(gca,'fontsize',20);
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
title('Eigenfunctions \phi (t^\prime)','fontsize',14,'fontweight','b');

%% construct X(t)
s = data;
x_n = zeros(N,1);
X = zeros(T,1);

for i = 1:N
   x_n(i) = dot(s, psi(i,:))/sig(i);
   X = X + x_n(i)* phi(i,:)';
end 

%% normalization
fac =T/sqrt(sum(X.^2));
X = X*fac;

%% receive signal r(t)
r = K*X;

r = r/fac;
% estimate/division
hat_r = round(r);
for i = 1:T2
   if hat_r(i)<-3
      hat_r(i) = -3;
   elseif hat_r(i) >3
      hat_r(i) = 3;
   end 
end 

e = s((T-T2):end) - hat_r((T-T2):end);
ne = length(find(e~=0));
BER = ne/T2

figure;
subplot(2,1,1)
plot([1:T],X);
grid on;
xlabel('t');
set(gca,'fontsize',20)
title('Transmit signal X(t)');
subplot(2,1,2)
plot([1:T],r(1:end),'.-');
hold on;
plot([1:T],s(1:end),'--');
grid on;
axis([1 100 -5 5])
xlabel('t');
legend('r(t)','s(t)');
set(gca,'fontsize',20)
title('Data s(t) and received signal r(t)');

%% test projection
% p = zeros(T,1);
% for i = 1:N
%    p = p+dot(s, psi(i,:))*(psi(i,:)');
% end
% s((100-T2):end) - p((100-T2):end);

%% Space interference

% 100 transmit antennas and 50 usrs, each usr has two antenna. 

S1 = 100; %100 users'antennas 
S2 = 100; %100 transmit antennas 
data2 = -3 + 2*randi([0,3],S1,1); %data serie

% generate space channel kernel.

Ks = zeros(S1, S2);

for i = 1:S1
   for j = 1:S2
      if i == j
         Ks(i,j) = 1;
      end 
      Ks(i,j) = randn();
   end
end 

figure;
mesh(Ks');
xlabel('u');
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
ylabel('u^\prime ');
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
% title('2D Kernel K(s1,s2)','fontsize',14,'fontweight','b');
set(gca,'fontsize',25)
%% 2D Generalized Mercer's theorem
[U2, Sig2, V2] = svd(Ks);
A2 = diag(Sig2);
sig2 = [];
for i = 1:length(A2)
   if A2(i) < 0.01 
      N2 = i-1; % N most contributed eigenfunctions
      break;
   end
   N2 = i;
   sig2(i) = A2(i);
end 

phi2 = zeros(N2,S2);
psi2 = zeros(N2,S1);

for i = 1:N2
   psi2(i,:) = U2(:,i);
   phi2(i,:) = V2(:,i);
end

% figure;
% subplot(2,1,1);
% [x3,y3] = meshgrid(1:size(psi2,2), 1:size(psi2,1));  % T is your table
% plot3(x3,y3,psi2);
% view(-8,67);
% xlabel('N','fontsize',12);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
% ylabel('u1','fontsize',12);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
% title('Eigenfunctions \psi (u1)','fontsize',14,'fontweight','b');
% subplot(2,1,2);
% [x4,y4] = meshgrid(1:size(phi2,2), 1:size(phi2,1));  % T is your table
% plot3(x4,y4,phi2);
% view(-8,67);
% xlabel('N','fontsize',12);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [-1.8 0 0.25])
% ylabel('u2','fontsize',12);
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [1.8 0 0.25])
% title('Eigenfunctions \phi (u2)','fontsize',14,'fontweight','b');

%% construct X(t)
s2 = data2;
x2 = zeros(N2,1);
X2 = zeros(S2,1);

for i = 1:N2
   x2(i) = dot(s2, psi2(i,:))/sig2(i);
   X2 = X2 + x2(i)* phi2(i,:)';
end 

%% receive signal r(t)
r2 = Ks*X2;

% estimate/division
hat_r2 = round(r2);
for i = 1:S1
   if hat_r2(i)<-3
      hat_r2(i) = -3;
   elseif hat_r2(i) >3
      hat_r2(i) = 3;
   end 
end 

e2 = s2(1:end) - hat_r2(1:end);
ne2 = length(find(e2~=0));
BER = ne2/S2

figure;
subplot(2,1,1)
stairs([1:S2],X2);
grid on;
xlabel('u');
title('Transmit signal X(u)');
set(gca,'fontsize',20)
subplot(2,1,2)
stairs([1:100],r2(1:end),'.-');
hold on;
stairs([1:100],s2(1:end),'--');
grid on;
axis([1 100 -5 5])
xlabel('u');
legend('r(u)','s(u)');
title('Data s(u) and received signal r(u)');
set(gca,'fontsize',20)

