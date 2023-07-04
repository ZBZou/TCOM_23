%%channel V2X
clear all
close all
%%
Nu = 10;
Nt = 10;
Bw = 20e6;
Ns = 64;
fs = Bw/Ns;
dt = 64/Bw;
N_snap = 4000;
fc = 5e9;
%% set tx
s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 5e9;                             % 2.4 GHz center frequency
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.sample_density = 10;                                 % Minimum possible sample density
% s.use_3GPP_baseline = 1;
l = qd_layout(s);
l.no_tx = 1; 
l.tx_position(:,:) = [ 0 ; 0 ; 10 ]; 
l.tx_array = qd_arrayant('3gpp-3d');
l.tx_array.no_elements = Nt;
l.tx_array.center_frequency = l.simpar.center_frequency;
l.tx_name = {'BS'};
l.tx_array.rotate_pattern(15,'z');
%% set rx
l.rx_array = qd_arrayant('vehicular'); 
l.rx_array.center_frequency = l.simpar.center_frequency;
l.no_rx = Nu;

%% track t1
% Sl = 'WINNER_UMa_C2_LOS';
Sn = '3GPP_38.901_UMa_NLOS';
l1 = 10;                       
o1 = 3*pi/4;        
v1 = 120;
t = l1/v1;
% t1.interpolate('time',10e-3,[],[],1);              
% t1.movement_profile = [ 0,0 ; 3,60 ; 6,150 ; 8,get_length(t1)]';
% dist1 = t1.interpolate_movement( 1e-3 );  
% t1.interpolate('time',10e-3,[],[],1);
% t1.calc_orientation;      


% Nu1 = randi(Nu);
% if Nu1 == Nu
%    Nu1 = Nu - 1;
% elseif Nu1 == 0
%    Nu1 = 1;
% end
display('set tracking');
for i = 1:Nu
   l.rx_array(1,i).no_elements = 1;
   l.rx_track(1,i) = qd_track.generate('street', l1, o1, 10, 30, 10, 50, 0);  
   len = l.rx_track(1,i).get_length;
   l.rx_track(1,i).set_speed(v1 + (-5+10*rand()));
%    l.rx_track(1,i).movement_profile = [1/10*t,1/10*len; 1/4*t, 1/2*len; 1/3*t, 1/2*len; t,len]';
   l.rx_track(1,i).name = append('UE',num2str(i));
   l.rx_track(1,i).initial_position = [200;0;1.5] + [10;10;0].*randn(3,1);
   l.rx_track(1,i).set_scenario(Sn,1, 2, 5,0);
%    l.rx_track(1,i).interpolate_positions( s.samples_per_meter ); 
%    l.rx_track(1,i).segment_index = [1, 2, 3];
%    l.rx_track(1,i).scenario = {Sn,Sn,Sn};
   l.rx_track(1,i).interpolate('time', dt);  
   l.rx_track(1,i).calc_orientation;
   l.rx_track(1,i).correct_overlap;                             % Adjust state change position
%    l.rx_track(1,i).no_snapshots = N_snap;
end
display('set done');

%% t2 

% l2 = 250;                   
% o2 = 3*pi/4;
% p2 = [ -200 ; 100 ; 2 ];        
% v2 = 120/3.6;
% 
% for i = Nu1+1:Nu
%    l.rx_array(1,i).no_elements = 1;
%    l.rx_track(1,i) = qd_track('linear',l2,o2);                                  % Set the rx-track
%    l.rx_track(1,i).set_speed(v2 + (-5+10*randn()));
%    l.rx_track(1,i).name = append('UE',num2str(i));
%    l.rx_track(1,i).initial_position = [0;0;1.5] + [2;5;0.1].*randn(3,1);
%    l.rx_track(1,i).set_scenario(Sn);
%    l.rx_track(1,i).interpolate_positions( s.samples_per_meter );  
%    l.rx_track(1,i).correct_overlap;                             % Adjust state change position
% end

%%
% No_seg = 8;
% Snaps_per_seg_1 = round(t1.no_snapshots/10);
% Snaps_per_seg_2 = round(t2.no_snapshots/10);
% 
% t1.segment_index = [0:1:No_seg-1]*Snaps_per_seg_1+1;
% t2.segment_index = [0:1:No_seg-1]*Snaps_per_seg_2+1;
% 
% S1 = '3GPP_38.901_UMa_NLOS';
% S2 = '3GPP_38.901_UMa_LOS';
% S3 = '3GPP_38.901_UMa_NLOS';
% S4 = '3GPP_38.901_UMa_LOS';
% S5 = '3GPP_38.901_UMa_LOS';
% S6 = '3GPP_38.901_UMa_NLOS';
% S7 = '3GPP_38.901_UMa_LOS';
% S8 = '3GPP_38.901_UMa_LOS';
% 
% 
% t1.scenario = {S1, S2, S3, S4, S5, S6, S7, S8};
% t1.interpolate_positions( s.samples_per_meter );
% 
% t2.scenario = {S1, S2, S3, S4, S5, S6, S7, S8};
% t2.interpolate_positions( s.samples_per_meter );
%% rx_rack
% 
% l.rx_track(1,1) = t1;
% l.rx_track(1,2) = t2;
% l.rx_track.correct_overlap;

%% para
p = l.init_builder;
gen_parameters(p);

%% channel generator
c = l.get_channels;                                    % Generate channel coefficients
cn = c.merge;
%% minimal snap
Nsnap = cn(1).no_snap; 
for i = 2:Nu
    Nsnap = [Nsnap, cn(i).no_snap];    
end
Minsnap = min(Nsnap);
%% Delay profile
% cn_quan = cn(1).quantize_delays(1/Bw);
% hc =  cn_quan.coeff; 
% hd = cn_quan.delay; 
% for i = 2:Nu
%    cn_quan_i = cn(i).quantize_delays(1/Bw);
%    h_i = cn_quan_i.coeff;
%    hc =  cat(1,hc(:,:,:,1:Minsnap),h_i(:,:,:,1:Minsnap));
%    hd_temp = cn_quan_i.delay;
%    hd = cat(1,hd(:,:,:,1:Minsnap),hd_temp(:,:,:,1:Minsnap));
% end
% Max_delay = max(hd,[],'all');
% NumDealyBins = Max_delay*Bw;
% h = zeros(Nu,Nt,NumDealyBins, Ns);
% for i = 1:Nu
%     for j = 1:Nu
%         for tau = 1:size(hc,3)
%             for t = 1:size(hc,4)
%                 bin_id = round(hd(i,j,tau,t)*Bw+1);
%                 h(i,j,bin_id,t) = hc(i,j,tau,t);
%             end
%         end
%     end
% end

%% Results
hf =  cn(1).fr(Bw,Ns); 
hf = hf(:,:,:,1:Minsnap);
for i = 2:Nu
   h_i = cn(i).fr(Bw,Ns);
   hf =  cat(1,hf(:,:,:,1:Minsnap),h_i(:,:,:,1:Minsnap));    
end
%%
h = ifft(hf,[],3);
%% plt ht
save_ht = 0;
ht_name = 'example_ht_v2x';
typet = 'time-delay';
plt_h_t = squeeze(h(1,1,:,1:Minsnap));
plot_h(plt_h_t, dt*1e3, 1e6/Bw, typet, save_ht, ht_name);
% figure;
% mesh([1:size(plt_h_t,1)]*1e6/Bw,[1:size(plt_h_t,2)]*dt, abs(plt_h_t)');
% xlabel('Delay(us)');
% ylabel('Time (s)');
% zlabel('Gain');
% view(-30,30);
% set(gca,'fontsize',14)
% title('Time-varying impulse reponse of Tx antenna 1 to user 1')
%% plt hf
save_ht = 0;
ht_name = 'example_hf_v2x';
typef = 'time-frequency';
plt_h_f = squeeze(hf(1,1,:,1:20));
plt_h_f = plt_h_f;
f_cor = (fc-Bw/2+[1:size(plt_h_f,1)]*Bw/size(plt_h_f,1))/1e6;
plot_h(plt_h_f, dt*1e3, f_cor, typef, save_ht, ht_name);
%% plt hs
plt_h_s = squeeze(h(:,:,1,1));
figure;
mesh(abs(plt_h_s)');
xlabel('Tx');
ylabel('Rx');
zlabel('Gain');
view(-30,30);
title('Spatial profile at time 1')
%% plt hst
save_ht = 0;
ht_name = 'example_ht_v2x';
typet = 'pdp';
plt_pdp = abs(squeeze(h(1,1,:,1:Minsnap))).^2;

%%
savefile = 1;

if savefile == 1
    save h_d10_v120.mat h -v7.3
end

%%
if savefile == 1
    save par_d10_v120.mat p Bw Ns fs dt Nu Nt fc v1 l1
end

%% power map
[ map,x_coords,y_coords] = l.power_map( '3GPP_38.901_UMa_LOS','quick',2,-100,300,-100,200,1.5 );
P = 10*log10( sum(sum(cat(3,map{:}),3),4));  
l.visualize([],[],0,1);                                   % Show BS and MT positions on the map
hold on
imagesc( x_coords, y_coords, P );                       % Plot the received power
hold off
axis([ -100 300 -100 200]);
caxis(max(P(:))+[-20 0]);
colmap = colormap;
colormap(colmap*0.5 + 0.5);
set(gca , 'layer' , 'top');

%% plot
% figure;
% plot (time , pow);
% xlabel('Time[s]');
% ylabel('Power[dB]');
% title('Path Gains');

% p1 = pow(1:100*r);
% p2 = pow(100*r+1 : 150*r);
% p3 = pow(150*r+1 : 240*r);
% p4 = pow(240*r+1 : 300*r);
% p5 = pow(300*r+1 : 400*r);
% p6 = pow(400*r+1 : end);
% 
% stat_p = [p1;p2;p3;p4;p5;p6];
% g  = [1*ones(length(p1),1);2*ones(length(p2),1);3*ones(length(p3),1);...
%    4*ones(length(p4),1);5*ones(length(p5),1);6*ones(length(p6),1)];
% % path_gains(:,:) = c.coeff(1,1,:,:);
% 
% figure
% boxplot(stat_p, g);
% xlabel('Segments');
% ylabel('Averaging path gains[dB]');
% title('Statistics');
% % figure
% % boxplot(10*log10(abs(path_gains))');
% 
% set(0,'defaultTextFontSize', 18)                      	% Default Font Size
% set(0,'defaultAxesFontSize', 18)                     	% Default Font Size
% set(0,'defaultAxesFontName','Times')               	    % Default Font Type
% set(0,'defaultTextFontName','Times')                 	% Default Font Type
% set(0,'defaultFigurePaperPositionMode','auto')       	% Default Plot position
% set(0,'DefaultFigurePaperType','<custom>')             	% Default Paper Type
% set(0,'DefaultFigurePaperSize',[14.5 7.8])            	% Default Paper Size
% 
% [ map,x_coords,y_coords] = l.power_map( '3GPP_38.901_UMa_LOS','quick',5,-500,500,-500,500,1.5 );
% P = 10*log10( sum(sum(cat(3,map{:}),3),4));                    % Total received power
% 
l.visualize([],[],0,1);                                   % Show BS and MT positions on the map
% hold on
% imagesc( x_coords, y_coords, P );                       % Plot the received power
% hold off
% axis([-300 300 -300 300])                               % Plot size
% caxis( max(P(:)) + [-20 0] )                            % Color range 
% colmap = colormap;
% colormap( colmap*0.5 + 0.5 );                           % Adjust colors to be "lighter"
% set(gca,'layer','top')                                  % Show grid on top of the map
p(1).visualize_clusters;