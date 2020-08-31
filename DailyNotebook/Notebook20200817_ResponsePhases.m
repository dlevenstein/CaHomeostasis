

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200817';

%% basic rate step (TTX)

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();

parms.kf_0 = 1e-4;      %(default: 1e-4) %CamKII-independent (i.e. baseline) phosphorylation rate
parms.kd_0 = 1e-7;  %(default: 1e-5)
parms.k_CaN = 0.11; %(default: 0.15)
parms.Ca_beta = -7.0; %(default: -7.05)
parms.s_beta = -35; %(default: -35)
parms.k_CamK = 0.15; %default: 0.1

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 150;
R_post = 25;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_WT')


% TTX Step: Increased Ca_PSP (TS)
% parms.Ca_PSP = 2;
% Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')
%parms.Ca_PSP = 1; %(Back to default)
%% Early phase: tau_b and kf_o
timeparms.maxT = 25*60;
taus_b = 140:20:500;
 kfs_0 = 0.5e-4:0.05e-4:2e-4;

%Reduced parms for tests
%  taus_b = 150:50:500;
%  kfs_0 = 0.5e-4:0.25e-4:2e-4;

clear simresults_delayphase
for bb = 1:length(taus_b)
    for kk = 1:length(kfs_0)
        bz_Counter(kk+(bb-1)*length(kfs_0),length(taus_b).*length(kfs_0),'sim')
    parms.tau_beta = taus_b(bb);
    parms.kf_0 = kfs_0(kk);
    simresults_delayphase(bb,kk) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
    end
end
 parms.tau_beta = 300; %Back to default
 parms.kf_0 = 1e-4;

%%
for bb = 1:length(taus_b)
    for kk = 1:length(kfs_0)
        counttime = simresults_delayphase(bb,kk).t_hr>0;
        try
        [peakheight(bb,kk),peakdelay(bb,kk)] = findpeaks(simresults_delayphase(bb,kk).m(counttime),'NPeaks',1,'SortStr','descend');
        catch
            %peakheight(bb,kk) = nan;
            %peakdelay(bb,kk) = nan;
            kk
            continue
        end
        temp =simresults_delayphase(bb,kk).t_hr(counttime);
        peakdelay(bb,kk) = temp(peakdelay(bb,kk));
    end
end
%%
figure
Plot_SynHomeo(simresults_delayphase(bb,kk),'FixYRange',1.75,'figwidth',2,'fignum',2,'title','test')

%%
figure
subplot(2,2,1)
imagesc(taus_b,kfs_0,peakheight')
colorbar
xlabel('Tau b');ylabel('k_f_0')
caxis([0 0.05])
axis xy
title('Peak Height')
subplot(2,2,2)
imagesc(taus_b,kfs_0,peakdelay')
xlabel('Tau b');ylabel('k_f_0')
colorbar
axis xy
title('Delay')

NiceSave('EarlyPhase',figfolder,[]) 


%% Peak Magnitude: k_cam
kms = [0.1:0.01:0.5];
%kms = [0.1:0.05:0.5];
timeparms.maxT = 36*60;

clear simresults
for bb = 1:length(kms)
    bz_Counter(bb,length(kms),'sim')
    parms.k_CamK = kms(bb);
    simresults_km(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
end
 parms.k_CamK = 0.15; %Back to original
 
%%
manip.blockB = @(t) 0.*(t>=t_TTX) + 0.*(t<t_TTX);
simresults_bblock = Run_SynHomeo(manip,parms,timeparms,'showfig',true);
manip = rmfield(manip,'blockB');
%%
simresults_kmall = bz_CollapseStruct(simresults_km,1);

%%
figure
subplot(2,2,1)
imagesc(simresults_km(bb).t_hr,kms,simresults_kmall.m)
axis xy
xlim([0 36])
xlabel('t (hr)');ylabel('k_C_a_M')

Plot_SynHomeo(simresults_bblock,'FixYRange',1.75,'figwidth',2,'fignum',2,'title','beta block')

NiceSave('CamKPeak',figfolder,[]) 

%%

