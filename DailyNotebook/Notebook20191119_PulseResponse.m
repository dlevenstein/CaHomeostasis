figfolder  = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Figures';


%% Parameters
kf_0 = 1;       %CamKII-independent (i.e. baseline) phosphorylation rate
kd_0 = 0.5;       %CaN-independent (i.e. baseline) dephosphorylation rate
k_CamK = 20;     %Maximal CamKII-mediated phoshorpylation rate
k_CaN = 30;      %Maximal CaN-mediated dephosphorylation rate (Free)

Ca_0 = -7.2;       %GluA1-independent (i.e. baseline) Calcium concentration
                    %(-6.854 at eqfrom SS)
                    
Ca_PSP = 0.1;     %Calcium from 1Hz PSPs with GluA1 fully phosphorylated

Ca_Kalpha = -5.5;   %Ca midpoint for CamKII alpha (data)
%Ca_Kbeta = 0.5;  %Ca midpoint for CamKII beta
Ca_Kdelta = 1;  %How much CamKIIbeta lowers the "threshold" of CamK
                %10x more sensitive (data)
Ca_CaN = -6.2;     %Ca midpoint for CaN (data)
Ca_beta = -6.8;    %Ca midpoint for CamKII alpha-beta transition

s_CamK = 8;    %Steepness of CamKII activation (data)
s_CaN = 4;      %Steepness of CaN activation (data)
%s_beta = -20;     %Steepness of CamKII alpha-beta transition
s_beta = -6;     %(Free)

tau_CamK = 1;   %Timescale of CamKII activation (data)
tau_CaN = 40;   %Timescale of CaN activation (data)
tau_beta = 300; %Timescale of CamKII alpha-beta transition

Ca_Kbeta = Ca_Kalpha-Ca_Kdelta;

%% 
f_upper = kf_0 + k_CamK.*0.8;
d_upper = kd_0 + k_CaN.*0.4;

f_lower = kf_0 + k_CamK.*0;
d_lower = kd_0 + k_CaN.*0.1;

A_upper = f_upper./(f_upper+d_upper);
A_lower = f_lower./(f_lower+d_lower);

Ruse = 20;
Ca_lower = Ca_0 + Ruse.*Ca_PSP.*A_lower;
Ca_upper = Ca_0 + Ruse.*Ca_PSP.*A_upper;
%%
camcolor = 'k';
cancolor = 'r';
betacolor = 'b';
Ca_X = linspace(Ca_0,-5,100);
R_X = linspace(0,50,100);
A_X = linspace(0,1,100);
[A_XY,R_XY] = meshgrid(A_X,R_X);

[M_XY,N_XY] = meshgrid(A_X,A_X);


figure

subplot(2,2,1)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_Kalpha,s_CamK),'color',camcolor,'linewidth',2)
    plot(Ca_X,Sigmoid(Ca_X,Ca_CaN,s_CaN),'color',cancolor,'linewidth',2)
    %plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)

    plot(Ca_X,Sigmoid(Ca_X,Ca_Kbeta,s_CamK),'--','color',camcolor,'linewidth',2)
    arrow('Start',[Ca_Kalpha-0.1,0.5],'Stop',[Ca_Kbeta+0.1,0.5],...
        'color',betacolor,'linewidth',1,'LineStyle',':','Length',10)
    legend('CamKII','CaN','\alpha -> \beta','location','northoutside')
axis tight
ylim([0 1])
    xlabel('Ca');ylabel('Activation')



subplot(4,2,5)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)
axis tight
ylim([0 1])
    xlabel('Ca');ylabel('% Beta')
    
    
subplot(3,3,3)
    imagesc(A_X,R_X,Ca_0 + R_XY.*Ca_PSP.*A_XY)
    ColorbarWithAxis([-7 -5],'logCa')
    axis xy
    xlabel('A (pGluA1)');ylabel('Rate');
    title('Calcium Entry')
    
subplot(3,3,6)
    Ainf = (kf_0 + k_CamK.*M_XY)./(kf_0 + k_CamK.*M_XY + kd_0 + k_CaN.*N_XY);
    
    imagesc(A_X,A_X,Ainf)
    ColorbarWithAxis([0 1],'A')
    axis xy
    xlabel('m (CaMKII)');ylabel('n (CaN)');
    title('GluA1 Phos')
    
subplot(3,3,9)
    tauA = 1./(kf_0 + k_CamK.*M_XY + kd_0 + k_CaN.*N_XY);
    imagesc(A_X,A_X,log10(tauA))
        colorbar

    ColorbarWithAxis([-2 0],'tau_A')
    LogScale('c',10)
    axis xy
    xlabel('m (CaMKII)');ylabel('n (CaN)');
    %title('GluA1 Phos')
NiceSave('ActivationFunctions',figfolder,[],'includeDate',true)

%%
maxT = 10*60;
dt = 0.01;
timesteps = -1000:dt:maxT;

%%
R_pre = 40;
R_pulse = [5 20];
t_app = [0 200];
pulsedur = 5;

R=R_pre.*ones(size(timesteps));

for pp = 1:length(R_pulse)
R(timesteps>t_app(pp) & timesteps<(t_app(pp)+pulsedur)) = R_pre+R_pulse(pp);
end


%% Updated parameterization
A = zeros(size(timesteps));
kf = zeros(size(timesteps));
kd = zeros(size(timesteps));
Ca = zeros(size(timesteps));
m = zeros(size(timesteps));
n = zeros(size(timesteps));
b = zeros(size(timesteps));

for tt = 2:length(timesteps)
    
    dAdt = kf(tt-1).*(1-A(tt-1)) - kd(tt-1).*A(tt-1);
    
    dmdt = (-m(tt-1) + Sigmoid(Ca(tt-1)+b(tt-1).*Ca_Kdelta,Ca_Kalpha,s_CamK))./tau_CamK;
    dndt = (-n(tt-1) + Sigmoid(Ca(tt-1),Ca_CaN,s_CaN))./tau_CaN;
    dbdt = (-b(tt-1) + Sigmoid(Ca(tt-1),Ca_beta,s_beta))./tau_beta;

    %dkfdt = (-kf(tt-1) + kf_0 + k_CamK.*Sigmoid(Ca(tt-1),(Ca_Kalpha.*(1-b(tt-1))+Ca_beta.*b(tt-1)),s_CamK))./tau_CamK;
    %dkddt = (-kd(tt-1) + kd_0 + k_CaN.*Sigmoid(Ca(tt-1),Ca_CaN,s_CaN))./tau_CaN;
    %dbetadt = (-b(tt-1) + Sigmoid(Ca(tt-1),Ca_beta,s_beta))./tau_beta;
    
    A(tt) = A(tt-1) + dAdt.*dt;
    m(tt) = m(tt-1) + dmdt.*dt;
    n(tt) = n(tt-1) + dndt.*dt;
    b(tt) = b(tt-1) + dbdt.*dt;
    
    Ca(tt) = Ca_0 + R(tt).*Ca_PSP.*A(tt);
    
    kf(tt) = kf_0 + k_CamK.*m(tt);
    kd(tt) = kd_0 + k_CaN.*n(tt);
    
end

%%
timwin = [-60 maxT]./60;
figure
subplot(5,1,1)
    plot(timesteps./60,R,'k')
    xlim(timwin)
    box off
    
    ylabel('Spike Rate')
subplot(5,1,2)
    plot(timesteps./60,A,'k')
    ylabel('pGluA1')
    xlim(timwin)
    ylim([0 1])
    box off
subplot(5,1,3)
    plot(timesteps./60,Ca,'k')
    ylabel('Ca')
    xlim(timwin)
    box off
    ylim([Ca_0 -2.5])
subplot(5,1,4)
    plot(timesteps./60,m,'r')
    hold on
    plot(timesteps./60,n,'b')
    ylabel('CamK/CaN Gates')
    legend('m','n')
    xlim(timwin)
    box off
    ylim([0 1])
subplot(5,1,5)
    plot(timesteps./60,b,'k')
    ylabel('% Beta')
    xlim(timwin)
    box off
    xlabel('T (hr)')
    ylim([0 1])

    
NiceSave('HebbianDynamics',figfolder,[],'includeDate',true)


%% Load the nullclines from XPP
nullfolder = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Notebook20191119/xppfigs';

nullnames = {'Baseline','Pulse','Recovery'};

for nn = 1:length(nullnames)
    [ null1{nn},null2{nn} ] = NullclinesFromXPP(fullfile(nullfolder,[nullnames{nn},'.dat'])); 
end

[ bifnlines ] = BifnFromXPP( fullfile(nullfolder,'nbbifn.dat') );

[ bifnlines_low ] = BifnFromXPP( fullfile(nullfolder,'nbbifn_lowrate.dat') );
%%
figure
for nn = 1:length(nullnames)
    subplot(3,3,nn)
    hold on
        plot(null1{nn}(:,1),null1{nn}(:,2),'r','linewidth',2)
        plot(null2{nn}(:,1),null2{nn}(:,2),'k','linewidth',2)
        box on
        xlabel('m');ylabel('A')
        title(nullnames{nn})
end

subplot(2,2,3)
hold on
plot(bifnlines(:,1),bifnlines(:,4),'k','linewidth',2)
plot(bifnlines(:,1),bifnlines(:,2),'k','linewidth',2)

plot(bifnlines(:,1),bifnlines(:,3),'k','linewidth',2)
plot(n(timesteps>0),b(timesteps>0))
xlim([0 1]);ylim([0 1])
box on
title('R = 40')
xlabel('n');ylabel('b')

subplot(2,2,4)
hold on
plot(bifnlines_low(:,1),bifnlines_low(:,2),'k','linewidth',2)
plot(bifnlines_low(:,1),bifnlines_low(:,3),'k','linewidth',2)
%plot(n(timesteps>0),b(timesteps>0))
xlim([0 1]);ylim([0 1])
box on
title('R = 20')

xlabel('n');ylabel('b')
NiceSave('PhasePlanes',figfolder,[],'includeDate',true)
