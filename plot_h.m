function plot_h(H, x_unit, y_unit, type, save_hst, filename)
if nargin < 2 
    time=1;
    delay=1;
    type = 'time-delay';
    save_hst = 0;
    filename = 'test';
elseif nargin < 4
    type = 'time-delay';
    save_hst = 0;
    filename = 'test';
elseif nargin < 5
    save_hst = 0;
    filename = 'test';    
end

if isscalar(x_unit)
    x_cor = [1:size(H,2)]*x_unit;
else
    x_cor = x_unit;
end 

if isscalar(y_unit)
    y_cor = [1:size(H,1)]*y_unit;
else
    y_cor = y_unit;
end 

if strcmp(type,'time-delay')
    figure;
    mesh(y_cor,x_cor, abs(H)');
    xlabel('Delay (\mus)');
    ylabel('Time (ms)');
    zlabel('Power');
    view(35,35);
    set(gca,'fontsize',16)
elseif strcmp(type,'time-frequency')
    figure;
    mesh(y_cor,x_cor, 10*log10(abs(H)'));
    xlabel('Frequency (Mhz)');
    ylabel('Time (ms)');
    zlabel('Power [dB]');
    view(35,35);
    set(gca,'fontsize',16)
elseif strcmp(type,'space-time')
    figure;
    mesh(y_cor,x_cor, (abs(H)'));
    xlabel('Users');
    ylabel('Delay (\mus)');
    zlabel('Power');
    view(135,35);
    set(gca,'fontsize',16)
elseif strcmp(type, 'space-space')
    figure;
    mesh(y_cor,x_cor, (abs(H)'));
    xlabel('Users');
    ylabel('Ants');
    zlabel('Power');
    view(35,35);
    set(gca,'fontsize',16)
elseif strcmp(type, 'kernel')
    figure;
    mesh(y_cor,x_cor, abs(H)');
    xlabel('u^\prime');
    ylabel('t^\prime (\mu s)');
    zlabel('Power');
    view(-135,35);
    xlim([1 size(H,1)]);
    set(gca,'fontsize',20);
    set(gca,'YTick',y_cor);
    set(gca,'YTickLabel',1000-y_cor);
end

[h, wd, ht] = tightfig();

if save_hst == 1
    name1 = append(filename, '.fig');
    name2 = append(filename, '.pdf');
    saveas(gcf, name1);
    exportgraphics(gcf, name2);
end
