figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210914';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.05e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 20;
%%
Ca = linspace(Ca_0,-5,100);


%%

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

%%

Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)
Ca_CamKII = -5.5;     %Ca midpoint for CaN (data)
s_CamKII = 8;      %Steepness of CaN activation (data)
Ca_b = -7.0;     %Ca midpoint for CaN (data)
s_b = -10;      %Steepness of CaN activation (data)
Ca_Kdelta = 1;
%Ca_Kdelta = 0;

kappa = 0.005;
kf_0 = 0.6e-3; 
k_CaN = kf_0./kappa;
k_CamK = 30.*k_CaN; %Free parameter

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
b = Sigmoid( Ca,Ca_b,s_b,1);
m = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII,s_CamKII,1);
kf = kf_0 + k_CamK.*m;
kd = k_CaN.*n;

A = kf./(kf+kd);
%R = 10.^(Ca - Ca_0)./(A.*Ca_PSP);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));

Ca_PSC0_TS = Ca_PSC0.*1.2;
mShift = 0.1;
m_TS = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII-mShift,s_CamKII,1);
kf_TS = kf_0 + k_CamK.*m_TS;
kd = k_CaN.*n;
A_TS = kf_TS./(kf_TS+kd);
R_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A_TS+1));

%k_CaNs
k_CaNs = k_CaN.*[0 0.02 0.06 0.2 0.6]';
A_kds = kf./(kf+k_CaNs.*n);
A_kds_TS = kf_TS./(kf_TS+k_CaNs.*n);
R_kd = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kds+1));
R_kd_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A_kds_TS+1));


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
    xlim([0 120])
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
    xlim([0 120])
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
    xlim([0 120])
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
    xlim([0 120])
    yrange = ylim(gca);    
    title('CaN-block')
    
NiceSave('RCa_CaNBlock',figfolder,[],'includeDate',true)


%% Simulate!

%% Block N: Post Only - Partial Block
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 120*60;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.3.*(t>=t_TTX) + 1.*(t<t_TTX);
%manip.Autophos = false;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
%parms.k_CamK = 0; %set to kCaN 6/22/21
parms.Ca_Kdelta = Ca_Kdelta;

kfs = [10.^-3.5 10.^-3 10.^-2.5 10.^-3.2];

clear simresults_NblockWT
for kk = 1:length(kfs)
    kk
    parms.kf_0 = kfs(kk);
    parms.k_CaN = parms.kf_0./kappa;
    parms.k_CamK = parms.k_CaN.*30;

    [simresults_NblockWT(kk),simparms] = Run_SynHomeo2(manip,parms,timeparms,...
        'showfig',true,'showfig_ActFun',true);
end
%%
manip_NMblock = manip;
manip_NMblock.blockM = @(t) 0.3.*(t>=t_TTX) + 1.*(t<t_TTX);
    [simresults_MNblock,simparms] = Run_SynHomeo2(manip_NMblock,parms,timeparms,...
        'showfig',true,'showfig_ActFun',true);
%%
figure
for kk = 1:length(kfs)
    Plot_SynHomeo(simresults_NblockWT(kk),'FixYRange',false,'fignum',2,'figwidth',2,'linecolor',[0.5 0.5 0.5])
    hold on
end
Plot_SynHomeo(simresults_NblockWT(kk),'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','k')
NiceSave('kf0',figfolder,[])
% Block N: TS
%%
parms_TS = parms;
parms_TS.Ca_Kalpha = simparms.Ca_Kalpha-mShift;     %TS mutation
parms_TS.Ca_PSP0 = Ca_PSC0_TS;
[simresults_NblockTS,simparms] = Run_SynHomeo2(manip,parms_TS,timeparms,'showfig',false);
%%
close all
Plot_SynHomeo(simresults_NblockTS,'manip',manip,'FixYRange',true,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_NblockWT(length(kfs)),'manip',manip,'FixYRange',true,'saveFig',figfolder,...
    'figname','BlockN_autophos_Partial','fignum',2,'figwidth',2,'title','80% Effective')

%%
%close all
%Plot_SynHomeo(simresults_NblockTS,'FixYRange',false,'fignum',2,'figwidth',2,'linecolor','r')
Plot_SynHomeo(simresults_MNblock,'manip',manip_NMblock,'FixYRange',true,'saveFig',figfolder,...
    'figname','BlockNM','fignum',2,'figwidth',2,'title','MNBlock')

%%
close all
Plot_SynHomeo(simresults_MNblock,'manip',manip_NMblock,'FixYRange',false,'saveFig',figfolder,...
    'fignum',2,'figwidth',2,'title','MNBlock')
Plot_SynHomeo(simresults_NblockWT(length(kfs)),'manip',manip,'FixYRange',false,'saveFig',figfolder,...
    'figname','BothExp','fignum',2,'figwidth',2,'linecolor','b')
