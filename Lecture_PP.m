% Copyright (C) Maximiliano Mollura and Riccardo Barbieri, 2019-2020.
% All Rights Reserved. See LICENSE.TXT for license details.
% maximiliano.mollura@polimi.it

%% Point Process Lecture 

clc, clear, close all

addpath(genpath(cd))

data_pth = fullfile(cd,'Original Signals');
ann_pth = fullfile(cd,'TiltData');

%% Import Data
lst = what(ann_pth);
lst = cell2mat(lst.mat);
% phase='Tilt6';
for i = 1%:size(lst,1)
    T = load(fullfile(ann_pth,lst(i,:)));
    tilt = T.(string(fieldnames(T)));
    load(fullfile(data_pth,lst(i,1:3)));
    
    ax = [];
    figure,
    ax = [ax,subplot(3,1,1)]; plot(Sig.tm,Sig.signal(:,1)), ylabel('RR [sec]')
    ax = [ax,subplot(3,1,2)]; plot(Sig.tm,Sig.signal(:,2)),...
    ax = [ax,subplot(3,1,3)]; plot(Sig.tm,Sig.signal(:,3)),...
    sgtitle(strcat('TIlt-Rest-StandUp Protocol - ',lst(i,:)))
    linkaxes(ax,'x')
    sgtitle(lst(i,:))
    
    ax=[];
    figure,
    ax = [ax,subplot(2,1,1)]; plot(tilt(:,1),tilt(:,2)), ylabel('RR [sec]')
    ax = [ax,subplot(2,1,2)]; plot(Sig.tm(Sig.tm>=tilt(1,1)&Sig.tm<=tilt(end,1)),...
        Sig.signal(Sig.tm>=tilt(1,1)&Sig.tm<=tilt(end,1),3)), ylabel('Tilt')
    linkaxes(ax,'x')
    sgtitle(lst(i,:))
end


%% Classical AR Modeling
ax = [];
R = tilt(100:200,1); % series of times of R-events [s]
Pmono = 6:13;
RRdep = diff(R)-mean(diff(R));
ttl = [];
% Order Selection
for i=1:length(Pmono)
mod=ar(RRdep,Pmono(i),'yw');
Ak(i)=aic(mod,'aic');
res=resid(RRdep,mod);
adt(i)=adtest(res.OutputData);
lil(i)=lillietest(res.OutputData);
ks(i)=kstest(res.OutputData);   
acf(i,:)=autocorr(res.OutputData);
ttl=[ttl;string(num2str(Pmono(i)))];
end
figure,
ax = [ax,subplot(1,3,1)]; plot(Pmono,Ak), title('Akaike'),ylabel('AIC'),xlabel('Order')
ax = [ax,subplot(1,3,2)]; plot(Pmono,adt), hold on, plot(Pmono,lil),
plot(Pmono,ks),legend('AD','Lil','KS'), title('ADTest'),ylabel('Null Hyp.'),xlabel('Order')
linkaxes(ax,'x')
subplot(1,3,3), plot(acf'),legend(ttl), ylabel('ACF'),xlabel('Lags'),title('ACF')

% Monovariate AR Modeling
Pmono = 8;
theta = aryule(RRdep,Pmono);
% theta = levinson(autocorr(RRdep),8);
[tot,f]=pyulear(RRdep,Pmono,512,1/mean(diff(R)));
% Spectral Analysis Mono
ax0=[];
fspect=figure;
ax0 = [ax0,subplot(1,2,1)]; plot(f,tot), title('Rest'),ylabel('PSD [sec^2/Hz]'),xlabel('Freqeuncy [Hz]')

R = tilt(400:500,1); % series of times of R-events [s]
RRdep = diff(R)-mean(diff(R));
Pmono = 6:13;
% Order
ttl=[];
for i=1:length(Pmono)
    mod=ar(RRdep,Pmono(i),'yw');
    Ak(i)=aic(mod,'aic');
    res=resid(RRdep,mod);
    adt(i)=adtest(res.OutputData);
    lil(i)=lillietest(res.OutputData);
    ks(i)=kstest(res.OutputData);
    acf(i,:)=autocorr(res.OutputData);
    ttl=[ttl;string(num2str(Pmono(i)))];
end
figure,
ax = [ax,subplot(1,3,1)]; plot(Pmono,Ak), title('Akaike'),ylabel('AIC'),xlabel('Order')
ax = [ax,subplot(1,3,2)]; plot(Pmono,adt), hold on, plot(Pmono,lil),
plot(Pmono,ks),legend('AD','Lil','KS'), title('ADTest'),ylabel('Null Hyp.'),xlabel('Order')
linkaxes(ax,'x')
subplot(1,3,3), plot(acf'),legend(ttl),ylabel('ACF'),xlabel('Lags'),title('ACF')

% Monovariate AR Modeling
Pmono=8;
theta=aryule(RRdep,Pmono);
[tot,f]=pyulear(RRdep,Pmono,512,1/mean(diff(R)));
% Spectral Analysis Mono
figure(fspect)
ax0 = [ax0,subplot(1,2,2)]; plot(f,tot), title('Tilt'),ylabel('PSD [sec^2/Hz]'),xlabel('Freqeuncy [Hz]')

linkaxes(ax0,'xy')


%% Batch Monovariate PP Modeling

R = tilt(100:200,1); % series of times of R-events [s]
Pmono = 8; %[4,6,8,10,12]

% Monovariate AR Modeling
[Thetap,Kappa,opt] = regr_likel(R, 'P', Pmono,'hasTheta0',1);
% Spectral Analysis Mono
meanRR = mean(diff(R)); % average RR interval
Var = meanRR / Kappa; % variance of an inverse Gaussian
[tot,~,~,f] = spectral(Thetap, Var, 1/meanRR);

ax = [];
f1 = figure;
ax = [ax,subplot(1,2,1)]; plot(f,tot), title('Rest'),ylabel('PSD [sec^2/Hz]'),xlabel('Freqeuncy [Hz]')

R = tilt(400:500,1); % series of times of R-events [s]
% Monovariave AR Modeling
[Thetap,Kappa,opt] = regr_likel(R, 'P', Pmono,'hasTheta0',1);
% Spectral Analysis Mono
meanRR = mean(diff(R)); % average RR interval
Var = meanRR / Kappa; % variance of an inverse Gaussian
[tot,~,~,f] = spectral(Thetap, Var, 1/meanRR);

figure(f1)
ax = [ax,subplot(1,2,2)];plot(f,tot), title('Tilt'),xlabel('Freqeuncy [Hz]')

linkaxes(ax,'xy')

%% Time-Variant Monovariate PP Modeling
R = tilt(:,1); % series of times of R-events [s]

Pmono = 9; %[4,6,8,10,12]
smpl = length(R); %[300,2000,length(R)]
[Mono.Thetap,Mono.Mu,Mono.Kappa,Mono.L,Mono.opt] = pplikel(R(1:smpl),'P',Pmono,'hasTheta0',1);
Mono.t = Mono.opt.t0 + (0:length(Mono.Mu)-1) * Mono.opt.delta;
Mono.Var = Mono.opt.meanRR.^3 ./ Mono.Kappa; % variance of an inverse Gaussian

ax = [];
figure;
ax = [ax,subplot(2,1,1)];
plot(R(2:end), diff(R))
hold on,plot(Mono.t, Mono.Mu), legend('RR','Mu')
ylabel('Tachogram [s]')
ax = [ax,subplot(2,1,2)];
plot(R(2:end), diff(R))
hold on,plot(Mono.t,Mono.opt.meanRR), legend('RR','MeanRR')
xlabel('time [s]')
ylabel('Tachogram [s]')

%% Hazard-Rate / KS-Plot / ACF
figure;
plot(Mono.t, Mono.L)
xlabel('t [s]'); ylabel('Lambda (hazard-rate) function');
[KSdistance,Z] = ks_plot(R(1:smpl), Mono.L, Mono.opt.delta);
[corr,threshold] = check_corr(Z);

%% Look at the PDFs
t=0:0.0005:1.3;
t=t*1000;
pd_rest_IG=makedist('InverseGaussian','mu',Mono.Mu(40000)*1000,'lambda',Mono.Kappa(40000)*1000);
pd_tilt_IG=makedist('InverseGaussian','mu',Mono.Mu(80000)*1000,'lambda',Mono.Kappa(80000)*1000);
% pd_rest_G=makedist('Normal','mu',Mono.Mu(40000)*1000,'sigma',sqrt((Mono.Mu(40000)*1000)^3/(Mono.Kappa(40000)*1000)));
% pd_tilt_G=makedist('Normal','mu',Mono.Mu(80000)*1000,'sigma',sqrt((Mono.Mu(80000)*1000)^3/(Mono.Kappa(80000)*1000)));
figure,plot(t,pdf(pd_rest_IG,t))
hold on,
% plot(t,pdf(pd_rest_G,t))
plot(t,pdf(pd_tilt_IG,t))
% plot(t,pdf(pd_tilt_G,t))
% legend('IG-rest','G-rest','IG-tilt','G-tilt')
legend('IG-rest','IG-tilt')

%% RR Spectra
undsmpl = 1; %
Thetap = Mono.Thetap(:,1:undsmpl:end);
Var = Mono.Var(:,1:undsmpl:end);
meanRR = Mono.opt.meanRR(:,1:undsmpl:end);
t = Mono.t(:,1:undsmpl:end);
[Mono.powLF, Mono.powHF, Mono.bal,~,Mono.powVLF,Mono.powTot] = hrv_indices(Thetap, Var, 1./meanRR);

ax = [];
figure,
ax = [ax,subplot(5,1,1)]; plot(t,Mono.powLF),hold on,
plot(t,Mono.powHF), legend('powLF','powHF')
ylabel('Power [sec^2]')
ax = [ax,subplot(5,1,2)];plot(t,Mono.powLF./(Mono.powTot-Mono.powVLF)), hold on,
plot(t,Mono.powHF./(Mono.powTot-Mono.powVLF)), legend('powLFn','powHFn')
ylabel('Power [n.u.]'), ylim([0 2])
ax = [ax,subplot(5,1,3)]; plot(t,Mono.bal), ylim([0 10])
ylabel('LF/HF');
ax = [ax,subplot(5,1,4)];yyaxis left, plot(t,Mono.Mu(1:undsmpl:end)), yyaxis right, plot(t,Var),legend('Mu','Var')
ax = [ax,subplot(5,1,5)]; plot(Sig.tm(Sig.tm>=tilt(1,1)&Sig.tm<=tilt(end,1)),...
    Sig.signal(Sig.tm>=tilt(1,1)&Sig.tm<=tilt(end,1),3)), ylabel('Tilt')
% xlabel('t [s]');
linkaxes(ax,'x')

%% Extract Features With Transitions
clear, close all, clc
Pmono = 10;
drct_ann = fullfile(cd,'TiltData');
lst = what(drct_ann);
lst = cell2mat(lst.mat);
DS = [];
for i=1:size(lst)
    
    % Load annotations
    T = load(fullfile(drct_ann,lst(i,:)));
    tilt = T.(string(fieldnames(T)));
    
    
    [Mu,Thetap,Kappa,L,opt.meanRR,opt.LogLikel] = deal([]);
    
    rest_t = [tilt(round(end/2),1)-150,tilt(round(end/2),1)-30]; % 2 mins
    tilt_t = [tilt(round(end/2),1)+30,tilt(round(end/2),1)+150]; % 2 mins
    [~,rest_s] = min(abs(rest_t-tilt(:,1)));
    [~,tilt_s] = min(abs(tilt_t-tilt(:,1)));
    % PP Mono
    EKGR = tilt(:,1);
    [Thetap,Mu,Kappa,L,opt] = pplikel(EKGR,'P',Pmono,'hasTheta0',1,'W',60);
    Var = opt.meanRR.^3 ./ Kappa; % variance of an inverse Gaussian
    [pLF, pHF, bal,~,pVLF,pTOT] = hrv_indices(Thetap, Var, 1./opt.meanRR);
    
    PP_t = [EKGR(1):opt.delta:EKGR(end)];
    rest_pp = PP_t>=rest_t(1) & PP_t<=rest_t(2); % 2 mins
    tilt_pp = PP_t>=tilt_t(1) & PP_t<=tilt_t(2); % 2 mins
    trans_pp = PP_t>rest_t(2) & PP_t<tilt_t(1); % 1 min
    
    % Features
    Feat = Extract_Features_all(tilt(rest_s(1):rest_s(2),2),...
        Mu(rest_pp),Var(rest_pp),pVLF(rest_pp),pLF(rest_pp),pHF(rest_pp),...
        bal(rest_pp),pTOT(rest_pp),"rest");
    Feat =[Feat;...
        Extract_Features_all(tilt(tilt_s(1):tilt_s(2),2),...
        Mu(tilt_pp),Var(tilt_pp),pVLF(tilt_pp),pLF(tilt_pp),pHF(tilt_pp),...
        bal(tilt_pp),pTOT(tilt_pp),string(fieldnames(T)))];
    Feat =[Feat;...
        Extract_Features_all(tilt(rest_s(2):tilt_s(1),2),...
        Mu(trans_pp),Var(trans_pp),pVLF(trans_pp),pLF(trans_pp),pHF(trans_pp),...
        bal(trans_pp),pTOT(trans_pp),strcat("Trans_",string(fieldnames(T))))];
    
    Feat.p_id = repmat(string(lst(i,1:3)),size(Feat,1),1);
    Feat = movevars(Feat,'p_id','Before','AVNN');
    DS = [DS;Feat];
    
    
end
mkdir('Datasets\')
save('Datasets\ML_DS_3','DS')