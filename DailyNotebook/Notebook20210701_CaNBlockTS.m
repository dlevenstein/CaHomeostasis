figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210701';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.2e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 4;
%%
Ca = linspace(Ca_0,-6,100);


%%

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

%%

Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)

kappa = 0.01;
kf_0 = 1e-3; 
k_CaN = kf_0./kappa;

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
kf = kf_0;
kd = k_CaN.*n;

A = kf./(kf+kd);
%R = 10.^(Ca - Ca_0)./(A.*Ca_PSP);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));

Ca_PSC0_TS = Ca_PSC0.*1.5;
R_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A+1));

%k_CaNs
k_CaNs = k_CaN.*[0 0.01 0.03 0.1 0.3]';
A_kds = kf./(kf+k_CaNs.*n);
R_kd = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kds+1));
R_kd_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A_kds+1));


%%
figure
subplot(2,2,1)
    plot((R),Ca,'k')
    hold on
    plot((R_none),Ca,'k--')
    plot((R_full),Ca,'k--')
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 150])
    yrange = ylim(gca);

subplot(2,2,2)
    plot(n,Ca)
    hold on
    plot(A,Ca)
    ylim(yrange)
    legend('n(Ca)','A(Ca)')
    box off
    xlabel('n, A')
   

    
subplot(3,3,7)
    plot((R),Ca,'k')
    hold on
    plot((R_kd),Ca,'k:')
    plot((R_none),Ca,'--','color',[0.5 0.5 0.5])
    plot((R_full),Ca,'--','color',[0.5 0.5 0.5])
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 150])
    yrange = ylim(gca);    
    title('CaN-block')
    
subplot(3,3,8)
    plot((R),Ca,'k')
    hold on
    plot((R_TS),Ca,'r')
    plot((R_kd(4,:)),Ca,'k:')
    plot((R_kd_TS(4,:)),Ca,'r:')
   % plot((R_none),Ca,'--','color',[0.5 0.5 0.5])
   % plot((R_full),Ca,'--','color',[0.5 0.5 0.5])

    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 150])
    yrange = ylim(gca);    
    title('CaN-block')
    
    
subplot(3,3,9)
    plot((R),A,'k')
    hold on
    plot((R_TS),A,'r')
    plot((R_kd(4,:)),A_kds(4,:),'k:')
    plot((R_kd_TS(4,:)),A_kds(4,:),'r:')
    ylabel('A')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 150])
    yrange = ylim(gca);    
    title('CaN-block')
    
NiceSave('RCa_CaNBlock',figfolder,[],'includeDate',true)

 %% Simulate CaEq as function of Kfo (keeping kappa constant)
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 120*60;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = kf_0./kappa;
parms.k_CamK = 0;
parms.k_CamK = 3.*parms.k_CaN; %Free parameter %do this once with 3, then 50
parms.Ca_Kdelta = 1.5;

R_pre = 60;
manip.rate = @(t) R_pre .* ones(size(t));

t_FK506 = 0;
manip.blockN = @(t) 0.1.*(t>=t_FK506) + 1.*(t<t_FK506);

TSmod_kf0 = linspace(0.5,1.5,9);
TSmod_CaPSP = linspace(1,2,9);
WT = [find(TSmod_kf0==1) find(TSmod_CaPSP==1)];
TS = [find(TSmod_kf0==1) find(TSmod_CaPSP==1.5)];
%kf0s_sim = kf_0.*TSmod_kf0;

clear equilvals_Ca equilvals_A equilvals_b equilvals_m equilvals_n
clear FKvvals_Ca FKvvals_A FKvvals_b FKvvals_m FKvvals_n
for kk = 1:length(TSmod_kf0)
    bz_Counter(kk,length(TSmod_kf0),'sim')
    
    parms.kf_0 = kf_0.*TSmod_kf0(kk);
    
    parfor pp = 1:length(TSmod_CaPSP)
        tempparms = parms;
        tempparms.Ca_PSP0 = Ca_PSC0 .* TSmod_CaPSP(pp);
        
        [simresults_kf0(kk,pp)] = Run_SynHomeo2(manip,tempparms,timeparms);
        equiltime = find(simresults_kf0(kk,pp).t_sec<0,1,'last');
        equilvals_Ca(kk,pp) = simresults_kf0(kk,pp).Ca(equiltime);
        equilvals_A(kk,pp) = simresults_kf0(kk,pp).A(equiltime);
        equilvals_b(kk,pp) = simresults_kf0(kk,pp).b(equiltime);
        equilvals_m(kk,pp) = simresults_kf0(kk,pp).m(equiltime);
        equilvals_n(kk,pp) = simresults_kf0(kk,pp).n(equiltime);
        
        FKvvals_Ca(kk,pp) = simresults_kf0(kk,pp).Ca(end);
        FKvvals_A(kk,pp) = simresults_kf0(kk,pp).A(end);
        FKvvals_b(kk,pp) = simresults_kf0(kk,pp).b(end);
        FKvvals_m(kk,pp) = simresults_kf0(kk,pp).m(end);
        FKvvals_n(kk,pp) = simresults_kf0(kk,pp).n(end);
        close all
    end
end

%%

figure
subplot(3,3,1)
imagesc(TSmod_CaPSP,TSmod_kf0,equilvals_Ca)
axis xy
colorbar
title('Ca')
crameri('berlin','pivot',equilvals_Ca(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,2)
imagesc(TSmod_CaPSP,TSmod_kf0,equilvals_A)
axis xy
colorbar
title('A')
crameri('berlin','pivot',equilvals_A(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,4)
imagesc(TSmod_CaPSP,TSmod_kf0,FKvvals_Ca)
axis xy
colorbar
title('Ca')
crameri('berlin','pivot',FKvvals_Ca(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,5)
imagesc(TSmod_CaPSP,TSmod_kf0,FKvvals_A)
axis xy
colorbar
title('A')
crameri('berlin','pivot',FKvvals_A(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,3)
imagesc(TSmod_CaPSP,TSmod_kf0,(equilvals_m.*parms.k_CamK)./kf_0)
axis xy
colorbar
title('m-mediated phos (*kf0)')
crameri('berlin','pivot',equilvals_m(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

NiceSave('TSspace_FK506_lowkCam',figfolder,[])
%%

close all
Plot_SynHomeo(simresults_kf0(TS(1),TS(2)),'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_kf0(WT(1),WT(2)),'FixYRange',false,'saveFig',figfolder,...
    'figname','FK506_lowkCam','fignum',2,'figwidth',2,'title','90% Effective')


 %% Simulate CaEq as function of Kfo (keeping kappa constant)
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 120*60;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = kf_0./kappa;
parms.k_CamK = 0;
parms.k_CamK = 50.*parms.k_CaN; %Free parameter %do this once with 3, then 50
parms.Ca_Kdelta = 1.5;

R_pre = 60;
manip.rate = @(t) R_pre .* ones(size(t));

t_FK506 = 0;
manip.blockN = @(t) 0.1.*(t>=t_FK506) + 1.*(t<t_FK506);

%Specified above...
% TSmod_kf0 = linspace(0.5,1.5,5);
% TSmod_CaPSP = linspace(1,2,5);
% WT = [find(TSmod_kf0==1) find(TSmod_CaPSP==1)];
% TS = [find(TSmod_kf0==1) find(TSmod_CaPSP==1.5)];
% %kf0s_sim = kf_0.*TSmod_kf0;

clear equilvals_Ca equilvals_A equilvals_b equilvals_m equilvals_n
clear FKvvals_Ca FKvvals_A FKvvals_b FKvvals_m FKvvals_n
for kk = 1:length(TSmod_kf0)
    bz_Counter(kk,length(TSmod_kf0),'sim')
    
    parms.kf_0 = kf_0.*TSmod_kf0(kk);
    
    parfor pp = 1:length(TSmod_CaPSP)
        tempparms = parms;
        tempparms.Ca_PSP0 = Ca_PSC0 .* TSmod_CaPSP(pp);
        
        [simresults_kf0(kk,pp)] = Run_SynHomeo2(manip,tempparms,timeparms);
        equiltime = find(simresults_kf0(kk,pp).t_sec<0,1,'last');
        equilvals_Ca(kk,pp) = simresults_kf0(kk,pp).Ca(equiltime);
        equilvals_A(kk,pp) = simresults_kf0(kk,pp).A(equiltime);
        equilvals_b(kk,pp) = simresults_kf0(kk,pp).b(equiltime);
        equilvals_m(kk,pp) = simresults_kf0(kk,pp).m(equiltime);
        equilvals_n(kk,pp) = simresults_kf0(kk,pp).n(equiltime);
        
        FKvvals_Ca(kk,pp) = simresults_kf0(kk,pp).Ca(end);
        FKvvals_A(kk,pp) = simresults_kf0(kk,pp).A(end);
        FKvvals_b(kk,pp) = simresults_kf0(kk,pp).b(end);
        FKvvals_m(kk,pp) = simresults_kf0(kk,pp).m(end);
        FKvvals_n(kk,pp) = simresults_kf0(kk,pp).n(end);
        close all
    end
end

%%

figure
subplot(3,3,1)
imagesc(TSmod_CaPSP,TSmod_kf0,equilvals_Ca)
axis xy
colorbar
title('Ca')
crameri('berlin','pivot',equilvals_Ca(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,2)
imagesc(TSmod_CaPSP,TSmod_kf0,equilvals_A)
axis xy
colorbar
title('A')
crameri('berlin','pivot',equilvals_A(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,4)
imagesc(TSmod_CaPSP,TSmod_kf0,FKvvals_Ca)
axis xy
colorbar
title('Ca')
crameri('berlin','pivot',FKvvals_Ca(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,5)
imagesc(TSmod_CaPSP,TSmod_kf0,FKvvals_A)
axis xy
colorbar
title('A')
crameri('berlin','pivot',FKvvals_A(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

subplot(3,3,3)
imagesc(TSmod_CaPSP,TSmod_kf0,(equilvals_m.*parms.k_CamK)./kf_0)
axis xy
colorbar
title('m-mediated phos (*kf0)')
crameri('berlin','pivot',equilvals_m(WT(1),WT(2)))
ylabel('kf0 mod');xlabel('CaPSP mod')

NiceSave('TSspace_FK506_highkCam',figfolder,[])
%%

close all
Plot_SynHomeo(simresults_kf0(TS(1),TS(2)),'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_kf0(WT(1),WT(2)),'FixYRange',false,'saveFig',figfolder,...
    'figname','FK506_highkCam','fignum',2,'figwidth',2,'title','90% Effective')
