

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210518';

%% Block M/N
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
manip.blockM = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

%parms.kf_0 = 1e-4;

[simresults_MNblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,...
    'FixYRange',true)

[simresults_MNblockWT_nobase,simparms] = Run_SynHomeo(manip,parms,timeparms,...
    'FixYRange',true,'blockBase',false)

figure
Plot_SynHomeo(simresults_MNblockWT,'FixYRange',true,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_MNblockWT_nobase,'FixYRange',true,'saveFig',figfolder,...
    'figname','BlockMN_WT','fignum',1,'figwidth',2)
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

[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
% Block N: TS
parms.Ca_PSP = 1.5;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false);

figure
Plot_SynHomeo(simresults_NblockTS,'FixYRange',1.5,'fignum',1,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_NblockWT,'FixYRange',1.5,'saveFig',figfolder,...
    'figname','BlockN_autophos','fignum',1,'figwidth',2)

%% Block N: Post Only - Partial Block
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.2.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = true;

[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
% Block N: TS
parms.Ca_PSP = 1.5;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false);


Plot_SynHomeo(simresults_NblockTS,'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_NblockWT,'FixYRange',false,'saveFig',figfolder,...
    'figname','BlockN_autophos_Partial','fignum',2,'figwidth',2,'title','80% Effective')




%% %%

%% Block N: Post Only - Change CaPSP
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.1.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = true;

%parms.kf_0 = 1e-4;

parms.Ca_PSP = 0.1;
[simresults_NblockWT_lowCaPSP,simparms] = Run_SynHomeo(manip,parms,timeparms)


% Block N: TS
parms.Ca_PSP = 0.15;
[simresults_NblockTS_lowCaPSP,simparms] = Run_SynHomeo(manip,parms,timeparms);

figure
Plot_SynHomeo(simresults_NblockTS_lowCaPSP,'FixYRange',false,'fignum',1,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_NblockWT_lowCaPSP,'FixYRange',false,'saveFig',figfolder,...
    'figname','BlockN_autophos_lowCaPSP','fignum',1,'figwidth',2)










%% Block N: Pre

parms.Ca_pre = -7;
parms.mini_max = 80;
parms.Ca_PSP = 1;
[preresults_NblockWT,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_WT','PreActFun','bell')

parms.Ca_PSP = 1.5;

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
manip.blockN = @(t) 0.1.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = false;

%parms.kf_0 = 1e-4;

[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_WT',...
    'FixYRange',true)


% Block N: TS
parms.Ca_PSP = 1.5;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_TS',...
    'FixYRange',true);

%% Block N: Pre

parms.Ca_pre = -7;
parms.mini_max = 80;
parms.Ca_PSP = 1;
[preresults_NblockWT,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','Pre_BlockN_WT','PreActFun','bell')

parms.Ca_PSP = 1.5;

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