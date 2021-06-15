figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210614';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration
                         %(-6.854 at eqfrom SS)
Ca_PSC0 = 0.3e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 3;
%%
Ca = linspace(Ca_0,-6.5,100);


%%

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

alphas = [0.5 2];
R_alphas =(10.^Ca - 10.^Ca_0)./(Ca_PSC0*(alphas'*1+1));

Ca_0s = [-8.5 -7.5]';
R_Ca0s =(10.^Ca - 10.^Ca_0s)./(Ca_PSC0.*(alpha.*0+1));

Ca_PSCs = [0.25e-8 1e-8]';
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
%%

kf_0 = 1e-4;       %CamKII-independent (i.e. baseline) phosphorylation rate
Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)
k_CaN = 0.11;


kf_0 = 1e-4; 
k_CaN = 1e-2;

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
kf = kf_0;
kd = k_CaN.*n;

A = kf./(kf+kd);
%R = 10.^(Ca - Ca_0)./(A.*Ca_PSP);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));

%kf0s
kf0s = [1e-3 1e-5]';
A_kf0 = kf0s./(kf0s+kd);
R_kf0 = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A_kf0+1));

%k_CaNs
k_CaNs = [1e-1 1e-3]';
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
%%




%% Simulate steady state calcium as f'n of rate

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();



timeparms.maxT = 0;
timeparms.preT = 4000;

R_pres = logspace(0.5,4,22);

for rr = 1:length(R_pres)
    bz_Counter(rr,length(R_pres),'sim')
    
R_pre = R_pres(rr);
manip.rate = @(t) R_pre .* ones(size(t));

[simresults,parms] = Run_SynHomeo(manip,parms,timeparms,'startON',false);
equilvals.Ca(rr) = simresults.Ca(end);
equilvals.b(rr) = simresults.b(end);
equilvals.m(rr) = simresults.m(end);
equilvals.n(rr) = simresults.n(end);
[simresults,parms] = Run_SynHomeo(manip,parms,timeparms,'startON',true);
equilvals.Ca_on(rr) = simresults.Ca(end);
equilvals.b_on(rr) = simresults.b(end);
equilvals.m_on(rr) = simresults.m(end);
equilvals.n_on(rr) = simresults.n(end);
close all
end


%%
Ca = linspace(Ca_0,-4,100);
n = Sigmoid( Ca,parms.Ca_CaN,parms.s_CaN,1);
kf = parms.kf_0;
kd = parms.k_CaN.*n;

A = kf./(kf+kd);
R = 10.^(Ca - parms.Ca_0)./(A.*parms.Ca_PSP);

kd = parms.kd_0;
A = kf./(kf+kd);
R_none = 10.^(Ca - parms.Ca_0)./(A.*parms.Ca_PSP);

%%
figure
subplot(2,2,1)
plot(log10(R_none),Ca,'k--')
hold on
plot(log10(R),Ca,'k')
plot(log10(R_pres),equilvals.Ca,'ko')
plot(log10(R_pres),equilvals.Ca_on,'ko')


legend('No CaN/CamKII','Theory: CaN Only','Simulation: With CamKII','location','northwest')
xlabel('PSP Rate');ylabel('Equilibrium Ca (log)')
box off
axis tight
xlim(log10(R_pres([1 end])))
LogScale('x',10)


subplot(2,2,2)
plot(n,Ca)
hold on
plot(Sigmoid( Ca,parms.Ca_beta,parms.s_beta,1),Ca)


subplot(4,2,5)
    hold on
    plot(log10(R_pres),(equilvals.b),'bo')
    plot(log10(R_pres),(equilvals.b_on),'bo')

    plot(log10(R_pres),(equilvals.m),'ro')
    plot(log10(R_pres),(equilvals.m_on),'ro')

    xlabel('PSP Rate');ylabel('CamKII Gates')
    box off
    axis tight
    xlim(log10(R_pres([1 end])))
    LogScale('x',10)

subplot(4,2,6)
    hold on
    plot(log10(R_pres),(equilvals.m),'ro')
    plot(log10(R_pres),(equilvals.m_on),'ro')

    xlabel('PSP Rate');ylabel('CamKII Gates')
    box off
    axis tight
    xlim(log10(R_pres([1 end])))
    LogScale('x',10)
    ylim([0 0.01])

subplot(4,2,7)
hold on

plot(log10(R_pres),equilvals.n,'ko')
plot(log10(R_pres),equilvals.n_on,'ko')


xlabel('PSP Rate');ylabel('CaN Gate')
box off
axis tight
xlim(log10(R_pres([1 end])))
LogScale('x',10)


NiceSave('RCa',figfolder,[],'includeDate',true)



%%
