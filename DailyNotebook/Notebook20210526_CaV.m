

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210526';



%% Block N: Post Only - Partial Block
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 100*60;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.1.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.Autophos = true;

parms.Ca_PSP = 0.5;
parms.Ca_V = 0.1;
[simresults_NblockWT,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false,'saveFig_AF',figfolder);
% Block N: TS
parms.Ca_PSP = 0.5;
parms.Ca_V = 0.15;
[simresults_NblockTS,simparms] = Run_SynHomeo(manip,parms,timeparms,'showfig',false);

close all
Plot_SynHomeo(simresults_NblockTS,'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_NblockWT,'FixYRange',false,'saveFig',figfolder,...
    'figname','BlockN_autophos_Partial','fignum',2,'figwidth',2,'title','90% Effective')


%%
Ca_PSPs = logspace(-2,0,9);
Ca_Vs = logspace(-2,-0.5,7);

clear simresults_screen
for vv = 1:length(Ca_Vs) 
    bz_Counter(vv,length(Ca_Vs),'v iteration')
    parfor pp = 1:length(Ca_PSPs)
        loopparms = parms;
        loopparms.Ca_PSP = Ca_PSPs(pp);
        loopparms.Ca_V = Ca_Vs(vv);
        [temp(pp)] = Run_SynHomeo(manip,loopparms,timeparms,'showfig',false);
    end 
    simresults_screen(:,vv) = temp;
end

%%
initIDX = find(simresults_screen(1,1).t_hr==-0.1);
equibIDX = length(simresults_screen(1,1).t_hr);
%figure
clear initCa equibCa initA equibA
for pp = 1:length(Ca_PSPs)
    for vv = 1:length(Ca_Vs)
        initA(pp,vv) = simresults_screen(pp,vv).A(initIDX); 
        equibA(pp,vv) = simresults_screen(pp,vv).A(equibIDX); 
        
        initCa(pp,vv) = simresults_screen(pp,vv).Ca(initIDX); 
        equibCa(pp,vv) = simresults_screen(pp,vv).Ca(equibIDX); 
        
%         subplot(length(Ca_PSPs),length(Ca_Vs),vv+(pp-1).*length(Ca_Vs))
%             plot(simresults_screen(pp,vv).t_hr,simresults_screen(pp,vv).A,'k')
%             xlim([0 24])
    end 
end
%%
ptest = 1;
vtest = 5;
figure
plot(simresults_screen(ptest,vtest).t_hr,simresults_screen(ptest,vtest).A,'k')
xlim([0 24])
%%
figure
subplot(2,2,1)
    imagesc(log10(Ca_Vs),log10(Ca_PSPs),log10(initA))
    xlabel('Ca_V');ylabel('Ca_P_S_P')
    LogScale('xy',10)
    colorbar
    axis xy
    caxis([-3 0])
    LogScale('c',10)
    
    
subplot(2,2,2)
    imagesc(log10(Ca_Vs),log10(Ca_PSPs),log10(equibA))
    colorbar
    axis xy
    xlabel('Ca_V');ylabel('Ca_P_S_P')
    LogScale('xy',10)
    caxis([-3 0])
    LogScale('c',10)
    
subplot(2,2,3)
    imagesc(log10(Ca_Vs),log10(Ca_PSPs),(initCa))
    xlabel('Ca_V');ylabel('Ca_P_S_P')
    LogScale('xy',10)
    colorbar
    axis xy
    %caxis([-3 0])
    %LogScale('c',10)
    
    
subplot(2,2,4)
    imagesc(log10(Ca_Vs),log10(Ca_PSPs),(equibCa))
    colorbar
    axis xy
    xlabel('Ca_V');ylabel('Ca_P_S_P')
    LogScale('xy',10)
    %caxis([-3 0])
    %LogScale('c',10)
    
NiceSave('CaParms',figfolder,[],'includeDate',true)
%% Reparameterize with RCaPSP and CaV/CaPSP
Ca_PSPs = logspace(-2,0.5,6);
Ca_VPSPs = logspace(-2,0,5);

clear simresults_screen_reparm
for vv = 1:length(Ca_VPSPs) 
    bz_Counter(vv,length(Ca_VPSPs),'v iteration')
    parfor pp = 1:length(Ca_PSPs)
        loopparms = parms;
        loopparms.Ca_PSP = Ca_PSPs(pp);
        loopparms.Ca_V = Ca_PSPs(pp).*Ca_VPSPs(vv);
        [temp(pp)] = Run_SynHomeo(manip,loopparms,timeparms,'showfig',false);
    end 
    simresults_screen_reparm(:,vv) = temp;
end

%%
initIDX = find(simresults_screen_reparm(1,1).t_hr==-0.1);
equibIDX = length(simresults_screen_reparm(1,1).t_hr);
figure
clear initCa equibCa initA equibA
for pp = 1:length(Ca_PSPs)
    for vv = 1:length(Ca_VPSPs)
        initA(pp,vv) = simresults_screen_reparm(pp,vv).A(initIDX); 
        equibA(pp,vv) = simresults_screen_reparm(pp,vv).A(equibIDX); 
        
        initCa(pp,vv) = simresults_screen_reparm(pp,vv).Ca(initIDX); 
        equibCa(pp,vv) = simresults_screen_reparm(pp,vv).Ca(equibIDX); 
        
        subplot(length(Ca_PSPs),length(Ca_VPSPs),vv+(pp-1).*length(Ca_VPSPs))
            plot(simresults_screen_reparm(pp,vv).t_hr,simresults_screen_reparm(pp,vv).A,'k')
            xlim([0 24])
    end 
end

%%
figure
subplot(2,2,1)
    imagesc(log10(Ca_VPSPs),log10(Ca_PSPs.*R_post),log10(initA))
    xlabel('Ca_V ratio');ylabel('R * Ca_P_S_P')
    LogScale('xy',10)
    colorbar
    axis xy
    caxis([-2 0])
    LogScale('c',10)
    
    
subplot(2,2,2)
    imagesc(log10(Ca_VPSPs),log10(Ca_PSPs.*R_post),log10(equibA))
    colorbar
    axis xy
    xlabel('Ca_V ratio');ylabel('R * Ca_P_S_P')
    LogScale('xy',10)
    caxis([-2 0])
    LogScale('c',10)
    
subplot(2,2,3)
    imagesc(log10(Ca_VPSPs),log10(Ca_PSPs.*R_post),(initCa))
    xlabel('Ca_V ratio');ylabel('R * Ca_P_S_P')
    LogScale('xy',10)
    colorbar
    axis xy
    %caxis([-3 0])
    %LogScale('c',10)
    
    
subplot(2,2,4)
    imagesc(log10(Ca_VPSPs),log10(Ca_PSPs.*R_post),(equibCa))
    colorbar
    axis xy
    xlabel('Ca_V ratio');ylabel('R * Ca_P_S_P')
    LogScale('xy',10)
    %caxis([-3 0])
    %LogScale('c',10)
    
NiceSave('CaParms_reparm',figfolder,[],'includeDate',true)