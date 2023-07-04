clear all;
close all;
%% complexity with users
Lu = [1:20];
Lt = [2000];

O1 = Lu.^3*Lt.^3;
O2 = ((Lu+Lt)/2).^5 + Lu.^2*Lt.^2;
O3 = Lt.*(Lu.^7 + Lu.^3).*factorial(Lu); 

f = figure;
semilogy(Lu, (O1),'LineWidth',2);
hold on;
semilogy(Lu, (O2),'LineWidth',2);
hold on;
semilogy(Lu, (O3),'LineWidth',2);
legend("HOGMT with Lemma 3", "HOGMT with HOSVD", "DPC",'location','southeast');
xlabel("Number of users");
ylabel("Complexity")
xlim([1,20]);
set(gca,'fontsize',16)
[h, wd, ht] = tightfig();
filename = 'complexity_1';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f, name1);
exportgraphics(f, name2);
%% complexity with snapshots
Lu = [10];
Lt = [1:4000];

O1 = Lu.^3*Lt.^3;
O2 = ((Lu+Lt)/2).^5 + Lu.^2*Lt.^2;
O3 = Lt.*(Lu.^7 + Lu.^3).*factorial(Lu); 
f = figure;
semilogy(Lt, (O1),'LineWidth',2);
hold on;
semilogy(Lt, (O2),'LineWidth',2);
hold on;
semilogy(Lt, (O3),'LineWidth',2);
legend("HOGMT with Lemma 3", "HOGMT with HOSVD", "DPC",'location','southeast');
xlabel("Number of data symbols");
ylabel("Complexity");
xlim([1,4000]);
set(gca,'fontsize',16)
[h, wd, ht] = tightfig();
filename = 'complexity_2';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f, name1);
exportgraphics(f, name2);