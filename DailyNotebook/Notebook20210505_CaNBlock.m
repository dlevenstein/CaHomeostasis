

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210505';

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
parms.Ca_PSP = 2;
Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')
%parms.Ca_PSP = 1; %(Back to default)

%% Block N: Post Only
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = true;

%parms.kf_0 = 1e-4;

[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_WT',...
    'FixYRange',true)


% Block N: TS
parms.Ca_PSP = 2;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_TS',...
    'FixYRange',true);

%% Block N: Pre

parms.Ca_pre = -7;
parms.mini_max = 80;
parms.Ca_PSP = 1;
[preresults_NblockWT,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_WT','PreActFun','bell')

parms.Ca_PSP = 2;

[preresults_NblockTS,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_TS','PreActFun','bell')


%% Plot A (mini Amplitude) and frequency
figure
    subplot(2,2,1)
        plot(simresults_NblockTS.t_hr,simresults_NblockTS.A,'r')
        hold on
        plot(simresults_NblockWT.t_hr,simresults_NblockWT.A,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel({'POST-ONLY MODEL','Mini Amplitude (A)'})
    subplot(2,2,2)
        plot(simresults_NblockTS.t_hr,simresults_NblockTS.R-R_post,'r')
        hold on
        plot(simresults_NblockWT.t_hr,simresults_NblockWT.R-R_post,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel('Mini Frequency (A)')

    subplot(2,2,3)
        plot(preresults_NblockTS.t_hr,preresults_NblockTS.A,'r')
        hold on
        plot(preresults_NblockWT.t_hr,preresults_NblockWT.A,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel({'PRE/POST MODEL','Mini Amplitude (A)'})
    subplot(2,2,4)
        plot(preresults_NblockTS.t_hr,preresults_NblockTS.R-R_post,'r')
        hold on
        plot(preresults_NblockWT.t_hr,preresults_NblockWT.R-R_post,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel('Mini Frequency (A)')
        
        NiceSave('NBlock',figfolder,[])
        
        
%% Block N: Post Only
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = false;

%parms.kf_0 = 1e-4;

[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_WT',...
    'FixYRange',true)


% Block N: TS
parms.Ca_PSP = 2;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_TS',...
    'FixYRange',true);

%% Block N: Pre

parms.Ca_pre = -7;
parms.mini_max = 80;
parms.Ca_PSP = 1;
[preresults_NblockWT,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_WT','PreActFun','bell')

parms.Ca_PSP = 2;

[preresults_NblockTS,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_TS','PreActFun','bell')


%% Plot A (mini Amplitude) and frequency
figure
    subplot(2,2,1)
        plot(simresults_NblockTS.t_hr,simresults_NblockTS.A,'r')
        hold on
        plot(simresults_NblockWT.t_hr,simresults_NblockWT.A,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel({'POST-ONLY MODEL','Mini Amplitude (A)'})
    subplot(2,2,2)
        plot(simresults_NblockTS.t_hr,simresults_NblockTS.R-R_post,'r')
        hold on
        plot(simresults_NblockWT.t_hr,simresults_NblockWT.R-R_post,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel('Mini Frequency (A)')

    subplot(2,2,3)
        plot(preresults_NblockTS.t_hr,preresults_NblockTS.A,'r')
        hold on
        plot(preresults_NblockWT.t_hr,preresults_NblockWT.A,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel({'PRE/POST MODEL','Mini Amplitude (A)'})
    subplot(2,2,4)
        plot(preresults_NblockTS.t_hr,preresults_NblockTS.R-R_post,'r')
        hold on
        plot(preresults_NblockWT.t_hr,preresults_NblockWT.R-R_post,'k')
        xlim([0 24])
        xlabel('t(hr)');ylabel('Mini Frequency (A)')
       
        NiceSave('NBlock',figfolder,'BlockAutoPhos')