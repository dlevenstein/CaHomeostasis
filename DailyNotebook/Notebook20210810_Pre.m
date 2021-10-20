figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210810';

%%
Ca_0 = -8.5;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.05e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 10;
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
Ca_b = -7.05;     %Ca midpoint for CaN (data)
s_b = -35;      %Steepness of CaN activation (data)
Ca_Kdelta = 1;

kappa = 0.01;
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


Ca_PSC0_TS = Ca_PSC0.*1.1;
m_TS = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII-0.1,s_CamKII,1);
kf_TS = kf_0 + k_CamK.*m_TS;
kd = k_CaN.*n;
A_TS = kf_TS./(kf_TS+kd);
R_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A_TS+1));


%Pre
R_evoked = 200;
Ca_pre = -6.9;
Ca_pre = -7.2;
s_pre = 25;
bellwidth = 0.3;
% Ca_pre = Ca_b;
% s_pre = s_b;
mini_max = 60;
mini_min = 10;
%p = Sigmoid( Ca,Ca_pre,s_pre,1).*(mini_max-mini_min)+mini_min;
p = Sigmoid(Ca,Ca_pre+bellwidth,-s_pre).*Sigmoid(Ca,Ca_pre,s_pre)*(mini_max-mini_min)+mini_min;
%%
figure
subplot(2,2,1)
    plot((R),Ca,'k')
    hold on
    plot((R_TS),Ca,'r')
    plot((R_none),Ca,'k--')
    plot((R_full),Ca,'k--')
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 250])
    yrange = ylim(gca);

subplot(2,2,2)
    plot(n,Ca)
    hold on
    plot(A,Ca)
    ylim(yrange)
    legend('n(Ca)','A(Ca)')
    box off
    xlabel('n, A')
   

subplot(2,2,3)
    plot(Ca,(R),'k')
    hold on
    plot(Ca,R_TS,'r')
    plot(Ca,(R_none),'k--')
    plot(Ca,R_full,'k--')
    plot(Ca,p,'b')
    plot(Ca,p+R_evoked,'b')
    plot(Ca,R_evoked.*ones(size(Ca)),'b--')
    ylabel('R')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('logCa');
    box off
    axis tight
    axis tight
     xlim([-8 -6])
     ylim([0 250])
    yrange = ylim(gca);

    
NiceSave('RCa',figfolder,[],'includeDate',true)


%% Simulate!

%% TTX
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 120*60;

R_pre = R_evoked;
R_post = 0;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
%manip.Autophos = false;
parms.Ca_pre = Ca_pre;
parms.s_pre = s_pre;
parms.bellwidth = bellwidth;
%parms.tau_pre = 300;
%parms.tau_beta = 300;
parms.mini_max = mini_max;
parms.mini_min = mini_min;

%parms.mini_max = 10;
%parms.mini_min = 10;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = k_CaN;
parms.k_CamK = k_CamK; %Free parameter
parms.Ca_Kdelta = Ca_Kdelta;


    [simresults,simparms] = Run_SynHomeo2_Pre(manip,parms,timeparms,...
        'showfig',true,'showfig_ActFun',true,'PreActFun','bell');
    %%
    parms_TS = parms;
    parms_TS.Ca_Kalpha = simparms.Ca_Kalpha-0.1;     %TS mutation
    %parms_TS.Ca_beta = simparms.Ca_beta+0.1; 
    parms_TS.Ca_PSP0 = Ca_PSC0_TS;
    [simresults_TS,simparms_TS] = Run_SynHomeo2_Pre(manip,parms_TS,timeparms,...
        'showfig',true,'showfig_ActFun',true,'PreActFun','bell');

%%

close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_TS,'FixYRange',false,'linecolor','r','saveFig',figfolder,...
    'figname','TTX_TSWT','fignum',2,'figwidth',2)
