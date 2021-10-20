figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210826';

%%
Ca_0 = -8.5;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
%Known solutin
Ca_PSC0 = 0.1e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 10;

% Ca_PSC0 = 0.1e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
% alpha = 6;
%%
Ca = linspace(Ca_0,-5,100);

%% VdC

Rs = linspace(0,100,100);
kvf = 0.004;
kvd_bar = 0.06;
kvd_l = 20;
kvd = kvd_bar.*exp(-Rs./kvd_l);
tau_V = 1./(kvf + kvd);
v_inf = kvf./(kvf + kvd);

kvd_l_TS = 8;
kvd_TS = kvd_bar.*exp(-Rs./kvd_l_TS);
tau_V_TS = 1./(kvf + kvd_TS);
v_inf_TS = kvf./(kvf + kvd_TS);



Ca_CamKII = -5.5;     %Ca midpoint for CaN (data)
s_CamKII = 8;      %Steepness of CaN activation (data)
Ca_b = -7.0;     %Ca midpoint for alpha-beta
s_b = -8;      %Steepness of alpha-beta activation
Ca_Kdelta = 1;

b = Sigmoid( Ca,Ca_b,s_b,1);
m = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII,s_CamKII,1);
mv = m' * v_inf;

m_alpha = Sigmoid( Ca,Ca_CamKII,s_CamKII,1);
mv_alpha = m_alpha' * v_inf;
m_beta = Sigmoid( Ca+Ca_Kdelta,Ca_CamKII,s_CamKII,1);
mv_beta = m_beta' * v_inf;


%%
figure
subplot(3,3,1)
plot(Rs,kvd,'k')
hold on
plot(Rs,kvf.*ones(size(Rs)),'b')
plot(Rs,kvd_TS,'r')
xlabel('R')
ylabel('k_v (min^-^1)')
legend('k_v_d','k_v_f')
box off


subplot(3,3,2)
plot(Rs,v_inf,'k')
hold on
plot(Rs,v_inf_TS,'r')
xlabel('R')
ylabel('v_i_n_f')
box off

subplot(3,3,3)
plot(Rs,tau_V,'k')
hold on
plot(Rs,tau_V_TS,'r')
xlabel('R')
ylabel('tau_v (min)')
box off

subplot(3,3,4)
imagesc(Ca,Rs,log10(mv_alpha)')
axis xy
xlabel('Ca')
ylabel('R')
ColorbarWithAxis([-5 0],'mv')
LogScale('c',10)
title('Alpha')

subplot(3,3,5)
imagesc(Ca,Rs,log10(mv_beta)')
axis xy
xlabel('Ca')
ylabel('R')
ColorbarWithAxis([-5 0],'mv')
LogScale('c',10)
title('Beta')

subplot(3,3,6)
imagesc(Ca,Rs,log10(mv)')
axis xy
xlabel('Ca')
ylabel('R')
ColorbarWithAxis([-5 0],'mv')
LogScale('c',10)
title('Alpha & Beta')

NiceSave('CdV_localization',figfolder,[],'includeDate',true)

%%

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

%%

Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)


kappa = 0.01;
kf_0 = 0.6e-3; 
k_CaN = kf_0./kappa;
k_CamK = 50.*k_CaN; %Free parameter
%k_CamK = 60.*k_CaN; %Free parameter

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
kf = kf_0 + k_CamK.*m;
kd = k_CaN.*n;

A = kf./(kf+kd);
%R = 10.^(Ca - Ca_0)./(A.*Ca_PSP);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));


Ca_PSC0_TS = Ca_PSC0.*1.1;
R_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A+1));


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
    xlim([0 200])
    yrange = ylim(gca);
    
NiceSave('RCa',figfolder,[],'includeDate',true)


%% Simulate!

%% TTX
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 120*60;

R_pre = 100;
R_post = 10;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
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

parms.kvf = kvf;
parms.kvd_bar = kvd_bar;
parms.kvd_l = kvd_l;
%Add: new model vs old

[simresults,simparms] = Run_SynHomeo2_CdV(manip,parms,timeparms,...
    'showfig',true,'showfig_ActFun',true);


parms_TS = parms;
parms_TS.kvd_l = kvd_l_TS;     %TS mutation
%parms_TS.Ca_beta = simparms.Ca_beta+0.1; 
parms_TS.Ca_PSP0 = Ca_PSC0_TS;
[simresults_TS,simparms_TS] = Run_SynHomeo2_CdV(manip,parms_TS,timeparms,...
    'showfig',true,'showfig_ActFun',true);


%%
close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',1,'figwidth',1)
Plot_SynHomeo(simresults_TS,'FixYRange',false,'linecolor','r','saveFig',figfolder,...
    'figname','TTX_TSWT','fignum',1,'figwidth',1)

%% With Presynaptic Rate Variation

R_res = @(t,R_eq,A,tau,f,phi) R_eq-A.*exp(-t./tau).*cos((t+phi).*2*pi./(f));

R_eq = 25; %Equilibrium rate
A = 25;     %Immediate post-TTX rate
tau = 26.*60;   %Decay time constant
f = 26.*60;     %Presynaptic oscillation frequency
phi = 0.*60;

%manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.rate = @(t) R_res(t,R_eq,A,tau,f,phi).*(t>=t_TTX) + R_pre.*(t<t_TTX);

figure
timepoints = linspace(0,72*60,100);
plot(timepoints./60,manip.rate(timepoints),'k')
%%
[simresults,simparms] = Run_SynHomeo2_CdV(manip,parms,timeparms,...
    'showfig',true,'showfig_ActFun',true);
 [simresults_TS,simparms_TS] = Run_SynHomeo2_CdV(manip,parms_TS,timeparms,...
     'showfig',true,'showfig_ActFun',true);
%%
close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_TS,'FixYRange',false,'linecolor','r','saveFig',figfolder,...
    'figname','TTX_TSWT_PreOsc','fignum',2,'figwidth',2)
