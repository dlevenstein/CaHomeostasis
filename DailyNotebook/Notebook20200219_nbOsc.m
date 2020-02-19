figfolder  = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Notebook20191126';


%% Parameters
kf_0 = 1;       %CamKII-independent (i.e. baseline) phosphorylation rate
kd_0 = 0.5;       %CaN-independent (i.e. baseline) dephosphorylation rate
k_CamK = 20;     %Maximal CamKII-mediated phoshorpylation rate
k_CaN = 30;      %Maximal CaN-mediated dephosphorylation rate (Free)
%k_CaN = 35;      %Maximal CaN-mediated dephosphorylation rate (Free)

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
maxT = 15*60;
dt = 0.01;
timesteps = -1000:dt:maxT;

%%
R_pre = 40;
R_pulse = [5 20 20 20];
t_app = [-250 -200 60 8*60];
pulsedur = 5;
R_post = 10;
t_TTX = [0 10*60];

R=R_pre.*ones(size(timesteps));

R(timesteps>t_TTX(1) & timesteps<t_TTX(2)) = R_post;

for pp = 1:length(R_pulse)
R(timesteps>t_app(pp) & timesteps<(t_app(pp)+pulsedur)) = ...
    R(timesteps>t_app(pp) & timesteps<(t_app(pp)+pulsedur))+R_pulse(pp);
end

%% Run it
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
timwin = [-300 maxT]./60;
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

    
NiceSave('PulseThenTTX',figfolder,[],'includeDate',true)


%% Load the nullclines from XPP
nullfolder = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Notebook20191126';

nullnames = {'UpperBranch','LowerBranch'};

for nn = 1:length(nullnames)
    [ null1{nn},null2{nn} ] = NullclinesFromXPP(fullfile(nullfolder,[nullnames{nn},'.dat'])); 
end

[ nbifnlines ] = BifnFromXPP( fullfile(nullfolder,'nbif_40.dat') );
[ mbifnlines ] = BifnFromXPP( fullfile(nullfolder,'mbif_40.dat') );
Cabifnlines(:,1) = nbifnlines(:,1);
[Cabifnlines(:,3)] = CaFromMN(mbifnlines(:,3),nbifnlines(:,3),40,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines(:,4)] = CaFromMN(mbifnlines(:,4),nbifnlines(:,4),40,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines(:,5)] = CaFromMN(mbifnlines(:,5),nbifnlines(:,5),40,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);

[ nbifnlines_20 ] = BifnFromXPP( fullfile(nullfolder,'nbif_20.dat') );
[ mbifnlines_20 ] = BifnFromXPP( fullfile(nullfolder,'mbif_20.dat') );
Cabifnlines_20(:,1) = nbifnlines_20(:,1);
[Cabifnlines_20(:,3)] = CaFromMN(mbifnlines_20(:,3),nbifnlines_20(:,3),20,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_20(:,4)] = CaFromMN(mbifnlines_20(:,4),nbifnlines_20(:,4),20,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_20(:,5)] = CaFromMN(mbifnlines_20(:,5),nbifnlines_20(:,5),20,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);

[ nbifnlines_50 ] = BifnFromXPP( fullfile(nullfolder,'nbif_50.dat') );
[ mbifnlines_50 ] = BifnFromXPP( fullfile(nullfolder,'mbif_50.dat') );
Cabifnlines_50(:,1) = nbifnlines_50(:,1);
[Cabifnlines_50(:,3)] = CaFromMN(mbifnlines_50(:,3),nbifnlines_50(:,3),50,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_50(:,4)] = CaFromMN(mbifnlines_50(:,4),nbifnlines_50(:,4),50,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_50(:,5)] = CaFromMN(mbifnlines_50(:,5),nbifnlines_50(:,5),50,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);

[ nbifnlines_30 ] = BifnFromXPP( fullfile(nullfolder,'nbif_30.dat') );
[ mbifnlines_30 ] = BifnFromXPP( fullfile(nullfolder,'mbif_30.dat') );
Cabifnlines_30(:,1) = nbifnlines_30(:,1);
[Cabifnlines_30(:,3)] = CaFromMN(mbifnlines_30(:,3),nbifnlines_30(:,3),30,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_30(:,4)] = CaFromMN(mbifnlines_30(:,4),nbifnlines_30(:,4),30,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
[Cabifnlines_30(:,5)] = CaFromMN(mbifnlines_30(:,5),nbifnlines_30(:,5),30,Ca_0,Ca_PSP,kf_0,k_CamK,kd_0,k_CaN);
%%
figure
for nn = 1:2
    subplot(3,3,nn)
    hold on
        plot(null1{nn}(:,1),null1{nn}(:,2),'r.','linewidth',2)
        plot(null2{nn}(:,1),null2{nn}(:,2),'b.','linewidth',2)
        box on
        xlabel('m');ylabel('n')
        title(nullnames{nn})
        xlim([0 1]);ylim([0 1])
end
NiceSave('PhasePlanes',nullfolder,[],'includeDate',true)
%%
%Ca_beta = 
%s_beta = 

figure


subplot(2,2,1)
hold on
plot(Cabifnlines(:,1),Cabifnlines(:,3),'k')
plot(Cabifnlines(:,1),Cabifnlines(:,4),'k:')
plot(Cabifnlines(:,1),Cabifnlines(:,5),'k')


plot(Cabifnlines_20(:,1),Cabifnlines_20(:,3),'k')
plot(Cabifnlines_20(:,1),Cabifnlines_20(:,4),'k:')
plot(Cabifnlines_20(:,1),Cabifnlines_20(:,5),'k')

plot(Cabifnlines_50(:,1),Cabifnlines_50(:,3),'k')
plot(Cabifnlines_50(:,1),Cabifnlines_50(:,4),'k:')
plot(Cabifnlines_50(:,1),Cabifnlines_50(:,5),'k')

plot(Cabifnlines_30(:,1),Cabifnlines_30(:,3),'k')
plot(Cabifnlines_30(:,1),Cabifnlines_30(:,4),'k:')
plot(Cabifnlines_30(:,1),Cabifnlines_30(:,5),'k')
plot(Sigmoid(Ca_X,Ca_beta,s_beta),Ca_X)
axis tight
xlim([0 1])
xlabel('b');ylabel('Ca (LogM)')

NiceSave('CaB_byRate',figfolder,[],'includeDate',true)


%%
Ca_beta_osc = -6.5;
s_beta_osc = -20;

figure


subplot(2,2,1)
hold on

plot(Cabifnlines_50(:,1),Cabifnlines_50(:,3),'k')
plot(Cabifnlines_50(:,1),Cabifnlines_50(:,4),'k:')
plot(Cabifnlines_50(:,1),Cabifnlines_50(:,5),'k')

plot(Cabifnlines(:,1),Cabifnlines(:,3),'k')
plot(Cabifnlines(:,1),Cabifnlines(:,4),'k:')
plot(Cabifnlines(:,1),Cabifnlines(:,5),'k')
plot(Sigmoid(Ca_X,Ca_beta_osc,s_beta_osc),Ca_X)
axis tight
xlim([0 1])
xlabel('b');ylabel('Ca (LogM)')
eee
NiceSave('Osc_CaB_byRate',figfolder,[],'includeDate',true)

%%
Ca_beta = Ca_beta_osc;
s_beta = s_beta_osc; 

%%
maxT = 72*60;
dt = 0.01;
timesteps = -1000:dt:maxT;

%%
R_pre = 50;
R_pulse = [0 0 0 0];
t_app = [-250 -200 60 8*60];
pulsedur = 5;
R_post = 40;
t_TTX = [0 maxT];

R=R_pre.*ones(size(timesteps));

R(timesteps>t_TTX(1) & timesteps<t_TTX(2)) = R_post;

for pp = 1:length(R_pulse)
R(timesteps>t_app(pp) & timesteps<(t_app(pp)+pulsedur)) = ...
    R(timesteps>t_app(pp) & timesteps<(t_app(pp)+pulsedur))+R_pulse(pp);
end


%% Run it
A = zeros(size(timesteps));
kf = zeros(size(timesteps));
kd = zeros(size(timesteps));
Ca = zeros(size(timesteps));
m = zeros(size(timesteps));
n = zeros(size(timesteps));
b = zeros(size(timesteps));

nred = 0.5; %new parameter - n reduces m

for tt = 2:length(timesteps)
    
    dAdt = kf(tt-1).*(1-A(tt-1)) - kd(tt-1).*A(tt-1);
    
    dmdt = (-m(tt-1) + (1-nred.*n(tt-1)).*Sigmoid(Ca(tt-1)+b(tt-1).*Ca_Kdelta,Ca_Kalpha,s_CamK))./tau_CamK;
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
timwin = [-300 maxT]./60;
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

    
NiceSave('RStep',figfolder,[],'includeDate',true)