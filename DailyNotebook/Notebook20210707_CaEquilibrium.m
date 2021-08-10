figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210707';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.2e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 4;
%%
Ca = linspace(Ca_0,-5,100);


%%

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

alphas = [0.5 2];
R_alphas =(10.^Ca - 10.^Ca_0)./(Ca_PSC0*(alphas'*1+1));

Ca_0s = [-8.5 -7.5]';
R_Ca0s =(10.^Ca - 10.^Ca_0s)./(Ca_PSC0.*(alpha.*0+1));

Ca_PSCs = [0.1e-8 0.3e-8]';
R_CaPSCs =(10.^Ca - 10.^Ca_0)./(Ca_PSCs.*(alpha.*0+1));
%%
figure
subplot(2,2,1)
    %plot((R),Ca,'k')
    hold on
    plot((R_none),Ca,'k')
    plot((R_full),Ca,'k--')

    ylabel('logCa')

    xlabel('R (Hz)')
    box off
    axis tight
    xlim([0 60])
    
    
subplot(3,3,7)
    %plot((R),Ca,'k')
    hold on
    plot((R_none),Ca,'k')
    plot((R_full),Ca,'k--')
    plot((R_alphas),Ca,'k:')

    ylabel('logCa')

    xlabel('R (Hz)')
    box off
    axis tight
    xlim([0 60])
    title('alpha')
    
    
subplot(3,3,8)
    %plot((R),Ca,'k')
    hold on
    plot((R_none),Ca,'k')
    %plot((R_full),Ca,'k--')
    plot((R_Ca0s),Ca,'k:')

    ylabel('logCa')

    xlabel('R (Hz)')
    box off
    axis tight
    xlim([0 60])
     title('Ca0')
     
subplot(3,3,9)
    %plot((R),Ca,'k')
    hold on
    plot((R_none),Ca,'k')
    %plot((R_full),Ca,'k--')
    plot((R_CaPSCs),Ca,'k:')

    ylabel('logCa')

    xlabel('R (Hz)')
    box off
    axis tight
    xlim([0 60])
     title('CaPSP')
     
    
    

NiceSave('RCa',figfolder,[],'includeDate',true)
%% Effects of kf0 and kCan ... kappa

Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)

kappa = 0.01;
kf_0 = 0.6e-3; 
k_CaN = kf_0./kappa;

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
kf = kf_0;
kd = k_CaN.*n;

A = kf./(kf+kd);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));

%kf0s
kf0s = [1e-2 1e-4]';
A_kf0 = kf0s./(kf0s+kd);
R_kf0 = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kf0+1));

%k_CaNs
k_CaNs = [1 1e-2]';
A_kds = kf./(kf+k_CaNs.*n);
R_kd = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kds+1));

%R1/2
Rs = linspace(0,60,60);
Ca_Ahalf = log10(10.^Ca_0 + Ca_PSC0.*Rs + Ca_PSC0.*alpha.*0.5.*Rs);
k_ratio = Sigmoid( Ca_Ahalf,Ca_CaN,s_CaN,1);

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
    plot((R_kf0),Ca,'k:')
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
    title('k_f_0')
    
subplot(3,3,8)
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
    title('k_C_a_N')
    
    
subplot(3,3,9)
    plot(k_ratio,Rs,'k')
    xlabel('k_f_0 / k_C_a_N')
    ylabel('R_A_0_._5')
    box off
    
    
NiceSave('RCa_A',figfolder,[],'includeDate',true)



%% Adding back in the m-gate (alpha)

Ca_CamKII = -5.5;     %Ca midpoint for CaN (data)
s_CamKII = 8;      %Steepness of CaN activation (data)

k_CamK = 50.*k_CaN; %Free parameter
m_alpha = Sigmoid( Ca,Ca_CamKII,s_CamKII,1);
kf_malpha = kf_0 + k_CamK.*m_alpha;

Ca_Kdelta = 1.5;
m_beta = Sigmoid( Ca+Ca_Kdelta,Ca_CamKII,s_CamKII,1);
kf_mbeta = kf_0 + k_CamK.*m_beta;

%Equilibrium beta
Ca_b = -7.05;     %Ca midpoint for CaN (data)
s_b = -35;      %Steepness of CaN activation (data)
b = Sigmoid( Ca,Ca_b,s_b,1);
m_ab = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII,s_CamKII,1);
kf_mab = kf_0 + k_CamK.*m_ab;

A_malpha = kf_malpha./(kf_malpha+kd);
R_malpha = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_malpha+1));
A_mbeta = kf_mbeta./(kf_mbeta+kd);
R_mbeta = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_mbeta+1));
A_mab = kf_mab./(kf_mab+kd);
R_mab = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_mab+1));

%kf0s
% kf0s = [1e-2 1e-4]';
% A_kf0 = kf0s./(kf0s+kd);
% R_kf0 = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kf0+1));

%k_CaNs
% k_CaNs = [1 1e-2]';
% A_kds = kf./(kf+k_CaNs.*n);
% R_kd = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kds+1));


%%
figure
subplot(2,2,1)
    plot((R),Ca,'k')
    hold on
    plot((R_mab),Ca,'r')
    plot((R_malpha),Ca,'b')
   % plot((R_mbeta),Ca,'g')
    
    plot((R_none),Ca,'k--')
    plot((R_full),Ca,'k--')
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 100])
    yrange = ylim(gca);
    legend('CaN-only','Alpha<->Beta CamKII','Alpha CamKII','location','southeast')

subplot(2,2,2)
    plot(n,Ca)
    hold on
    plot(m_alpha,Ca)
    plot(m_beta,Ca)
    plot(A,Ca)
    plot(A_malpha,Ca)
    plot(A_mbeta,Ca)
    ylim(yrange)
    legend('n(Ca)','m_a(Ca)','m_a(Ca)','A(Ca)')
    box off
    xlabel('n, A')
   
subplot(3,3,7)
    plot((R),Ca,'k')
    hold on
    plot((R_kf0),Ca,'k:')
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
    title('k_f_0')
    
subplot(3,3,8)
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
    title('k_C_a_N')
    
    
subplot(3,3,9)
    plot(k_ratio,Rs,'k')
    xlabel('k_f_0 / k_C_a_N')
    ylabel('R_A_0_._5')
    box off
    
    
NiceSave('RCa_A_mgate',figfolder,[],'includeDate',true)



%% Simulate steady state calcium as f'n of rate

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 0;
timeparms.preT = 6000;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = 0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = kf_0./kappa;
parms.k_CamK = k_CamK;
parms.Ca_Kdelta = Ca_Kdelta;

R_pres = logspace(0,2,5);

for rr = 1:length(R_pres)
    bz_Counter(rr,length(R_pres),'sim')
    
    R_pre = R_pres(rr);
    manip.rate = @(t) R_pre .* ones(size(t));

    parms.Ca_A = 0;
    parms.kf_0 = kf_0;
    [simresults,parms] = Run_SynHomeo2(manip,parms,timeparms);
    equilvals.Ca(rr) = simresults.Ca(end);
    equilvals.b(rr) = simresults.b(end);
    equilvals.m(rr) = simresults.m(end);
    equilvals.n(rr) = simresults.n(end);
    close all

    parms.kf_0 = 100;
    parms.Ca_A = alpha.*Ca_PSC0;
    [simresults,parms] = Run_SynHomeo2(manip,parms,timeparms);
    equilvals_Full.Ca(rr) = simresults.Ca(end);
    equilvals_Full.b(rr) = simresults.b(end);
    equilvals_Full.m(rr) = simresults.m(end);
    equilvals_Full.n(rr) = simresults.n(end);
    close all

    parms.kf_0 = kf_0;
    parms.Ca_A = alpha.*Ca_PSC0;
    [simresults,parms] = Run_SynHomeo2(manip,parms,timeparms);
    equilvals_Eq.Ca(rr) = simresults.Ca(end);
    equilvals_Eq.b(rr) = simresults.b(end);
    equilvals_Eq.m(rr) = simresults.m(end);
    equilvals_Eq.n(rr) = simresults.n(end);
    close all
end

%%
figure
plot(simresults.t_hr,simresults.Ca,'k')
%%
figure
    subplot(2,2,1)
        plot((R),Ca,'k')
        hold on
        plot((R_none),Ca,'k--')
        plot((R_full),Ca,'k--')
        plot(R_pres,equilvals.Ca,'o')
        plot(R_pres,equilvals_Full.Ca,'o')
        plot(R_pres,equilvals_Eq.Ca,'o')
        ylabel('logCa')
        %plot(log10(R),A,'r')
        %LogScale('x',10)
        xlabel('R');
        box off
        axis tight
        axis tight
        xlim([0 60])
        yrange = ylim(gca);


 %% Simulate CaEq as function of Kfo (keeping kappa constant)
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 0;
timeparms.preT = 6000;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = kf_0./kappa;
parms.k_CamK = 0;

R_pre = 30;
manip.rate = @(t) R_pre .* ones(size(t));

kf0s_sim = [1e-2 1e-3 1e-4 1e-5];

for rr = 1:length(kf0s_sim)
    bz_Counter(rr,length(kf0s_sim),'sim')
    
    parms.kf_0 = kf0s_sim(rr);
    parms.k_CaN = kf0s_sim(rr)./kappa;
    [simresults_kf0,parms] = Run_SynHomeo2(manip,parms,timeparms);
    equilvals_kf0.Ca(rr) = simresults_kf0.Ca(end);
    equilvals_kf0.b(rr) = simresults_kf0.b(end);
    equilvals_kf0.m(rr) = simresults_kf0.m(end);
    equilvals_kf0.n(rr) = simresults_kf0.n(end);
    close all

end

%%
figure
subplot(3,3,7)
    plot((R),Ca,'k')
    hold on
    plot((R_kf0),Ca,'k:')
    plot((R_none),Ca,'--','color',[0.5 0.5 0.5])
    plot((R_full),Ca,'--','color',[0.5 0.5 0.5])
    plot(R_pre,equilvals_kf0.Ca,'o')
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 60])
    yrange = ylim(gca);
    title('k_f_0')