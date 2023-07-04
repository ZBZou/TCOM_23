close all
clear all

%% Para
Nu = 8;
Nt = 8;
Ns = 2000;
Bw = 20e6;
fc = 5e9;
Nf = 64;
df = Bw/Nf;
dt = 1/df;
ds = 2000;
visual_position = 0;
speed = 120/3.6;
save = 0;
%% Use the winner2.wimparset function to create a WINNER II model parameter set.

cfgwim = winner2.wimparset;
cfgwim.RandomSeed = 31; % Set the rng seed for repeatability
cfgwim.NumTimeSamples = ds;
cfgwim.CenterFrequency = fc;

%% Define antenna arrays for one BS and two MS.
Bs  = winner2.AntennaArray('UCA', Nt, 0.02);  
Ms = [];
for i = 1:Nu
    Ms = [Ms, winner2.AntennaArray('ULA', 1, 0.01)];
end

%% Create system layout by using the winner2.layoutparset function.
MSIdx = [2:Nu+1]; 
BSIdx = {1}; 
K = Nu; 
rndSeed = 5;
cfgLayout = winner2.layoutparset(MSIdx,BSIdx, ...
    K,[Bs,Ms],[],rndSeed);
%% Senario setting
cfgLayout.ScenarioVector = 12*ones(1,Nu);
cfgLayout.PropagConditionVector = zeros(1,Nu);
for i = 1:Nu
    cfgLayout.Stations(i+1).Velocity = [speed;speed;0].*randn(3,1);
end
%% Visualize BS and MS positions.
if visual_position == 1
    BSPos  = cfgLayout.Stations(cfgLayout.Pairing(1,1)).Pos;
    
    figure;
    hold on;
    grid on;
    plot3(BSPos(1),BSPos(2),BSPos(3),'bo');
    for i = 1:Nu
        MSPos  = cfgLayout.Stations(cfgLayout.Pairing(2,i)).Pos;
        plot3(MSPos(1),MSPos(2),MSPos(3),'ro');
    end
    
    xlim([0 500]);
    ylim([0 500]);
    zlim([0 35]);
    xlabel('X-position (m)');
    ylabel('Y-position (m)');
    zlabel('Elevation (m)');
    view(-30,30);
    legend('BS','UEs','Location','northeast');
end
%% Generate channel coefficients
[Coef, Path_delay, Cond] = winner2.wim(cfgwim,cfgLayout);
for i = 2: Ns/ds
    [Coef_i, Path_delay, Cond] = winner2.wim(cfgwim,cfgLayout,Cond);
    Coef = cellfun(@(x,y) cat(4,x,y),Coef,Coef_i,'UniformOutput',false);
end
%% Merge links
H =Coef{1};
for i = 2:Nu
   H_i = Coef{i};
   H =  cat(1,H,H_i);    
end
%% 
WINNERChan = comm.WINNER2Channel(cfgwim,cfgLayout);
chanInfo = info(WINNERChan);

%% Convert to H_t
Max_delay = max(Path_delay(:,end));
NumDealyBins = Max_delay/cfgwim.DelaySamplingInterval;
h = zeros(Nu,Nt,NumDealyBins, Ns);
for i = 1:Nu
    for j = 1:Nt
        for tau = 1:size(Path_delay,2)
            for t = 1:Ns
                bin_id = round(Path_delay(i,tau)/cfgwim.DelaySamplingInterval)+1;
                h(i,j,bin_id,t) = H(i,j,tau,t);
            end
        end
    end
end

%% Plot H_t
plt_h_t = squeeze(h(1,1,:,:));
figure;
mesh([1:size(plt_h_t,1)]*cfgwim.DelaySamplingInterval,[1:size(plt_h_t,2)]*Cond.delta_t(1), abs(plt_h_t)');
%%
if save == 0
    save winner_v2x.mat h -v7.3
end