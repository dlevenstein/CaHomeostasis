figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200831';

%%
Ca_0 = -7.5;       %GluA1-independent (i.e. baseline) Calcium concentration
                            %(-6.854 at eqfrom SS)
Ca_PSP = 1;     %Calcium per 1Hz PSPs with GluA1 fully phosphorylated
kf_0 = 1e-4;       %CamKII-independent (i.e. baseline) phosphorylation rate
Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)

k_CaN = 0.11;


%%
Ca = linspace(Ca_0,-6,100);


%%
n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
kf = kf_0;
kd = k_CaN.*n;

A = kf./(kf+kd);
R = 10.^(Ca - Ca_0)./(A.*Ca_PSP);

%%
figure
subplot(2,2,1)
plot(log10(R),Ca,'k')
LogScale('x',10)

subplot(2,2,2)
plot(n,Ca)

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
