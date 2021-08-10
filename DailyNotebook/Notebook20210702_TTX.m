figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210702';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.2e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 4;
%%
Ca = linspace(Ca_0,-6.5,100);


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
    xlim([0 60])
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
    xlim([0 60])
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
    xlim([0 60])
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
    xlim([0 60])
    yrange = ylim(gca);    
    title('CaN-block')
    
NiceSave('RCa_CaNBlock',figfolder,[],'includeDate',true)


%% Simulate!

%% TTX
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 120*60;

R_pre = 60;
R_post = 10;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
%manip.Autophos = false;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = 1e-3;
parms.k_CaN = parms.kf_0./kappa;
parms.k_CamK = 50.*parms.k_CaN; %Free parameter
parms.Ca_Kdelta = 1.5;


[simresults,simparms] = Run_SynHomeo2(manip,parms,timeparms,...
    'showfig',true,'showfig_ActFun',true);


parms_TS = parms;
Ca_PSC0_TS = Ca_PSC0.*1.5;
parms_TS.Ca_PSP0 = Ca_PSC0_TS;
[simresults_TS,simparms] = Run_SynHomeo2(manip,parms_TS,timeparms,...
    'showfig',true,'showfig_ActFun',true);

%%

close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_TS,'FixYRange',false,'linecolor','r','saveFig',figfolder,...
    'figname','TTX_TSWT','fignum',2,'figwidth',2)
