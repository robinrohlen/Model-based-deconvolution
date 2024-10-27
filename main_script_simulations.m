%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example script to simulate a force signal and estimate the twitch
% parameters using:
% 1) model-based deconvolution
% 2) spike-triggered averaging (STA)
%
% See preprint for more info: https://doi.org/10.1101/2024.05.14.594072
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all

T=60; % Total simulation time in seconds
fs=1e3; % Sample rate
t=linspace(0,T,T*fs); % Time vector
ramp=[5 50 5]; % up-stable-down in seconds
max_exc_lvl=2.5; % Percentage % of max exc

% Motoneuron params
n_mn=200; % Number of motoneurons in the pool
RR=30; % Range of recruitment threshold values
g_e=1; % gain of the excitatory drive-firing rate relationship
MFR=8; % Hz - Minimum firing rate
PFR1=35; % Hz - peak firing rate of the first motoneuron
PFRD=10; % Hz - desired difference in peak firing rates between the first and last units
ISICoV=0.15; % Inter-spike interval coefficient of variation (15%)

a=log(RR)./n_mn;
RTE=exp(a.*(1:n_mn)); % Recruitment threshold excitation

PFRi=PFR1-PFRD*(RTE./RTE(end)); % Peak firing rate of motoneuron i
E_max=RTE(end)+(PFRi(end)-MFR)./g_e; % maximum excitation

E_t=[linspace(0,(max_exc_lvl/100)*E_max,ramp(1)*fs) (max_exc_lvl/100)*E_max*ones(1,ramp(2)*fs) flip(linspace(0,(max_exc_lvl/100)*E_max,ramp(3)*fs))]; % Excitatory drive

for i=1:n_mn
    t_thresh=E_t-RTE(i); % Above this thresh - fire
    find_t_thresh=find(t_thresh>=0); % Which samples are associated with firing
    if isempty(find_t_thresh)
        continue;
    end
    curr_t=t(find_t_thresh(1));
    exc_diff=t_thresh(find_t_thresh(1)); % Time point for the first impulse
    tmp=max(1./(g_e*(exc_diff)+MFR),1./PFRi(i));
    t_imp{i}(1)=(ISICoV*tmp)*randn(1)+curr_t;
    iter=1; % Firing counter
    while curr_t+1/MFR<=t(find_t_thresh(end)) % while sample point is below last sample point associated with firing
        tmp=max(1./(g_e*(exc_diff)+MFR),1./PFRi(i));
        t_imp{i}(iter+1)=(ISICoV*tmp)*randn(1)+tmp+curr_t;
        iter=iter+1; % Next firing
        % Find next threshold after current time
        findThresholdSample=t(find_t_thresh)-t_imp{i}(iter);
        minInd=find(findThresholdSample>0);
        if isempty(minInd)
            continue
        else
            minInd=minInd(1);
        end
        exc_diff=t_thresh(find_t_thresh(minInd));
        % Set current time
        curr_t=t(find_t_thresh(minInd));
    end
end

%%% Generate twitch params
RP=100;

Tlead=randi([5 5],1,size(t_imp,2)); % [5 15] in paper, but cannot be reliably estimated
Ttw=randi([1000 1000],1,size(t_imp,2)); % fix length, uses tukey window for compact support

Tc=150*(1./exp((log(RP)./n_mn).*(1:n_mn))).^(1/2.861353); % 150 to 30 ms
Thr=250*(1./exp((log(RP)./n_mn).*(1:n_mn))).^(1/2.861353); % 250 to 50 ms
Fmax=exp((log(RP)./n_mn).*(1:n_mn));

Tc=Tc(1:size(t_imp,2));
Thr=Thr(1:size(t_imp,2));
Fmax=Fmax(1:size(t_imp,2));

%%% Generate force based on spike trains and twitch params

t=linspace(0,round(1e3*max(cell2mat(t_imp)))+max(Tlead)+max(Ttw),(fs*10^(-3))*(round(1e3*max(cell2mat(t_imp)))+max(Tlead)+max(Ttw)));
t=ceil(t);
Fsum=cell(size(t_imp,2));
force=zeros(size(t));

for MU=1:size(t_imp,2)
    Ftmp=zeros(size(t));
    for ind=1:length(t_imp{MU})
        % Generate twitch using model (using Timp=0 and Tlead=0)
        [~,F]=RaikovaForceTwitch5p(fs,round(1e3*t_imp{MU}(ind)),Tlead(MU),Tc(MU),Thr(MU),Ttw(MU),Fmax(MU));
        % Padd with zeros using Timp and Tlead size(t,2)
        Fpadded=zeros(size(t));
        [~,idx]=min(abs(t-(round(1e3*t_imp{MU}(ind))+Tlead(MU))));
        Fpadded(idx:(idx+length(F)-1))=F;
        % Linear superposition of twitches
        Ftmp=Ftmp+Fpadded;
    end
    Fsum{MU}=Ftmp;
    force=force+Ftmp;
end

%%% Add noise + high-pass filter
snr_const=20;

[b,a] = butter(3,4/(fs/2),'high'); % 4 Hz was found to be the best
noise_const=mad(force)*10^(-snr_const/20);
force_noise=force+noise_const.*randn(size(force));
force_filt=filtfilt(b,a,force_noise);

%%% Estimate params (deconv optim)

est_params_sta=cell(1,size(t_imp,2));
est_params_deconv=cell(1,size(t_imp,2));
true_params_deconv=cell(1,size(t_imp,2));

options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',@patternsearch,'Display','off');

iter=5;
nvars=3; % in the original paper 4
lb=[20 4 0.5]; % in the original paper [0 20 4 0.5];
ub=[180 6 150]; % in the original paper [20 180 6 150];

parfor MUind=1:size(t_imp,2)
    disp([num2str(MUind),'/',num2str(size(t_imp,2))])
    locs=round(1e3*t_imp{MUind});
    x=zeros(iter,nvars);

    fun=@(theta) norm(twitch_optim(theta,force_filt,locs,fs)-force_filt);

    for ind=1:iter
        [xtmp,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options);
        x(ind,:)=xtmp;
    end
    xmean=mean(x);

    win=0:200;
    twitch_STA=[];
    for spike_num=1:size(t_imp{MUind},2)
        twitch_STA=cat(1,twitch_STA,force_noise(round(1e3*t_imp{MUind}(spike_num))+win));
    end

    est_params_sta{MUind}=cat(1,est_params_sta{MUind},[find(mean(twitch_STA)-min(mean(twitch_STA))==max(mean(twitch_STA)-min(mean(twitch_STA)))) max(mean(twitch_STA)-min(mean(twitch_STA)))]);
    est_params_deconv{MUind}=cat(1,est_params_deconv{MUind},[0 xmean(1) xmean(1)*xmean(2)/3 xmean(3)]);
    true_params_deconv{MUind}=cat(1,true_params_deconv{MUind},[Tlead(MUind) Tc(MUind) Thr(MUind) Fmax(MUind)]);

end

save(['example',num2str(max_exc_lvl),'.mat'])

% Plots
tmp1=cell2mat(true_params_deconv);
tmp2=cell2mat(est_params_deconv);
tmp3=cell2mat(est_params_sta);

cmap=lines(3);

h=tiledlayout(2,2);set(gcf,'units','points','position',[315,113,936,729]);

nexttile;
hold on;
p1=plot(100*(tmp1(2:4:end)-tmp2(2:4:end))./tmp1(2:4:end),'o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12);
p2=plot(100*(tmp1(2:4:end)-tmp3(1:2:end))./tmp1(2:4:end),'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12);
p3=plot(100*(tmp1(2:4:end)-(median(tmp2(2:4:end))+tmp3(1:2:end)-median(tmp3(1:2:end))))./tmp1(2:4:end),'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerSize',12);
plot(xlim,[0 0],'k:')
hold off;
xlabel('Motor unit #');
ylabel('Ground truth - Estimated (%)')
title('Contraction time')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
ylim([-100 100]);
l=legend([p1 p2 p3],{'Deconv','STA','STA (bias adj.)'},'Orientation','horizontal');

nexttile;
hold on;
plot(100*(tmp1(4:4:end)-tmp2(4:4:end))./tmp1(4:4:end),'o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12);
plot(100*(tmp1(4:4:end)-tmp3(2:2:end))./tmp1(4:4:end),'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12);
plot(100*(tmp1(4:4:end)-(median(tmp2(4:4:end))+tmp3(2:2:end)-median(tmp3(2:2:end))))./tmp1(4:4:end),'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerSize',12);
plot(xlim,[0 0],'k:')
hold off;
xlabel('Motor unit #');
ylabel('Ground truth - Estimated (%)')
title('Twitch force')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
ylim([-100 100]);

nexttile;
hold on;
plot(tmp1(2:4:end)-tmp2(2:4:end),'o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12);
plot(tmp1(2:4:end)-tmp3(1:2:end),'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12);
plot(tmp1(2:4:end)-(median(tmp2(2:4:end))+tmp3(1:2:end)-median(tmp3(1:2:end))),'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerSize',12);
plot(xlim,[0 0],'k:')
hold off;
xlabel('Motor unit #');
ylabel('Ground truth - Estimated (ms)')
title('Contraction time')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);

nexttile;
hold on;
plot(tmp1(4:4:end)-tmp2(4:4:end),'o','Color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'MarkerSize',12);
plot(tmp1(4:4:end)-tmp3(2:2:end),'o','Color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'MarkerSize',12);
plot(tmp1(4:4:end)-(median(tmp2(4:4:end))+tmp3(2:2:end)-median(tmp3(2:2:end))),'o','Color',cmap(3,:),'MarkerFaceColor',cmap(3,:),'MarkerSize',12);
plot(xlim,[0 0],'k:')
hold off;
xlabel('Motor unit #');
ylabel('Ground truth - Estimated (ms)')
title('Twitch force')
set(gca,'TickDir','out');set(gcf,'color','w');set(gca,'FontSize',16);
