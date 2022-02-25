figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20220131';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
%Known solutin
Ca_PSC0 = 0.05e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 20;

% Ca_PSC0 = 0.1e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
% alpha = 6;
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
   

    

    
NiceSave('RCa',figfolder,[],'includeDate',true)


%% Simulate!

%% TTX
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 40*60;
timeparms.preT = 120*60;


R_pre = 100;
R_post = 10;
t_TTX = 0;

%R_baseline = 100;
R_pulse = 150;
t_pulse = [-600 1500]+t_TTX;
pulsedur = 5;

manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX) + ...
    R_pulse.*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse.*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur));

%manip.Autophos = false;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = k_CaN;
parms.k_CamK = k_CamK; %Free parameter
parms.Ca_beta = Ca_b;     %Ca midpoint for alpha-beta
parms.s_beta = s_b;      %Steepness of alpha-beta activation
parms.Ca_Kdelta = Ca_Kdelta;



[simresults,simparms] = Run_SynHomeo2(manip,parms,timeparms,...
    'showfig',true);

parms_BB = parms;
parms_BB.Ca_Kdelta = 0;     %TS mutation
%parms_TS.Ca_beta = simparms.Ca_beta+0.1; 
[simresults_BB,simparms_BB] = Run_SynHomeo2(manip,parms_BB,timeparms,...
    'showfig',true,'showfig_ActFun',true);



close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_BB,'FixYRange',false,'linecolor','b','saveFig',figfolder,...
    'figname','TTX_PulseResponse_BB','fignum',2,'figwidth',2)

%% With Presynaptic Rate Variation

% R_res = @(t,R_eq,A,tau,f,phi) R_eq-A.*exp(-t./tau).*cos((t+phi).*2*pi./f);
% 
% R_eq = 25; %Equilibrium rate
% A = 20;     %Immediate post-TTX rate
% tau = 50.*60;   %Decay time constant
% f = 35.*60;     %Presynaptic oscillation frequency
% phi = 5.*60;
% 
% 
% %manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
% manip.rate = @(t) R_res(t,R_eq,A,tau,f,phi).*(t>=t_TTX) + R_pre.*(t<t_TTX);
% 
% [simresults,simparms] = Run_SynHomeo2(manip,parms,timeparms,...
%     'showfig',true,'showfig_ActFun',true);
% [simresults_BB,simparms_BB] = Run_SynHomeo2(manip,parms_BB,timeparms,...
%     'showfig',true,'showfig_ActFun',true);
% %%
% close all
% Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
% Plot_SynHomeo(simresults_BB,'FixYRange',false,'linecolor','b','saveFig',figfolder,...
%     'figname','TTX_BB_PreOsc','fignum',2,'figwidth',2)

% Plot_SynHomeo(simresults,'FixYRange',false,'linecolor','k','saveFig',figfolder,...
%     'figname','TTX_WT_PreOsc','fignum',2,'figwidth',2)