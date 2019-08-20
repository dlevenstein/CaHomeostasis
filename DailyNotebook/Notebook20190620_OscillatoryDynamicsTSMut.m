figfolder  = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Figures';


%% Parameters
kf_0 = 1;       %CamKII-independent (i.e. baseline) phosphorylation rate
kd_0 = 0.5;       %CaN-independent (i.e. baseline) dephosphorylation rate
k_CamK = 10;     %Maximal CamKII-mediated phoshorpylation rate
k_CaN = 20;      %Maximal CaN-mediated dephosphorylation rate

Ca_0 = 0;       %GluA1-independent (i.e. baseline) Calcium concentration    
Ca_PSP = 0.7;     %Calcium from 1Hz PSPs with GluA1 fully phosphorylated
Ca_Kalpha = 3.5;   %Ca midpoint for CamKII alpha
%Ca_Kbeta = 0.5;  %Ca midpoint for CamKII beta
Ca_Kdelta = 3;  %How much CamKIIbeta lowers the "threshold" of CamK
Ca_CaN = 2.25;     %Ca midpoint for CaN
Ca_beta = 0.8;    %Ca midpoint for CamKII alpha-beta transition

s_CamK = 10;    %Steepness of CamKII activation
s_CaN = 1;      %Steepness of CaN activation
s_beta = -20;     %Steepness of CamKII alpha-beta transition

tau_CamK = 1;   %Timescale of CamKII activation
tau_CaN = 40;   %Timescale of CaN activation
tau_beta = 500; %Timescale of CamKII alpha-beta transition

Ca_Kbeta = Ca_Kalpha-Ca_Kdelta;

%% 
f_upper = kf_0 + k_CamK.*1;
d_upper = kd_0 + k_CaN.*0.3;

f_lower = kf_0 + k_CamK.*0;
d_lower = kd_0 + k_CaN.*0.2;

A_upper = f_upper./(f_upper+d_upper);
A_lower = f_lower./(f_lower+d_lower);

Ruse = 2.7;
Ca_lower = Ca_0 + Ruse.*Ca_PSP.*A_lower;
Ca_upper = Ca_0 + Ruse.*Ca_PSP.*A_upper;
%%
camcolor = 'k';
cancolor = 'r';
betacolor = 'b';
Ca_X = linspace(0,5,100);

figure

subplot(2,2,1)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_Kalpha,s_CamK),'color',camcolor,'linewidth',2)
    plot(Ca_X,Sigmoid(Ca_X,Ca_CaN,s_CaN),'color',cancolor,'linewidth',2)
    plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)

    plot(Ca_X,Sigmoid(Ca_X,Ca_Kbeta,s_CamK),'--','color',camcolor,'linewidth',2)
    arrow('Start',[Ca_Kalpha-0.1,0.5],'Stop',[Ca_Kbeta+0.1,0.5],...
        'color',betacolor,'linewidth',1,'LineStyle',':','Length',10)
    legend('CamKII','CaN','\alpha -> \beta','location','southeast')

    xlabel('Ca');ylabel('Activation')


subplot(2,2,2)
    hold on
    plot(Ca_X,kf_0 + k_CamK.*Sigmoid(Ca_X,Ca_Kalpha,s_CamK),'color',camcolor,'linewidth',2)
    plot(Ca_X,kd_0 + k_CaN.*Sigmoid(Ca_X,Ca_CaN,s_CaN),'color',cancolor,'linewidth',2)

    plot(Ca_X,kf_0 + k_CamK.*Sigmoid(Ca_X,Ca_Kbeta,s_CamK),'--','color',camcolor,'linewidth',2)
    arrow('Start',[Ca_Kalpha-0.1,0.5.*(k_CamK)+kf_0],...
        'Stop',[Ca_Kbeta+0.1,0.5.*(k_CamK)+kf_0],...
        'color',betacolor,'linewidth',1,'LineStyle',':','Length',10)
    %legend('CamKII','CaN','alpha/beta','location','southeast')
    xlabel('Ca');ylabel('De/phorphorylation rate')

subplot(4,2,6)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)

    xlabel('Ca');ylabel('% Beta')
    
subplot(4,2,8)
    hold on
    plot(Ca_X,Ca_Kdelta.*Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)
    axis tight
    xlabel('Ca');ylabel('Camred')

NiceSave('ActivationFunctions',figfolder,[],'includeDate',true)

%%
maxT = 7500;
dt = 0.1;
timesteps = 0:dt:maxT;

%%
R_pre = 5;
R_ttx = 2.7;
t_app = 2500;
R_spk = 0;%exp(randn(size(timesteps)));
R_spk = movmean(R_spk,1./dt);
%R_spk = 0;
R=R_pre.*ones(size(timesteps));
R(timesteps>t_app) = R_ttx;
R=R+R_spk;
A = zeros(size(timesteps));
kf = zeros(size(timesteps));
kd = zeros(size(timesteps));
beta = zeros(size(timesteps));
Ca = zeros(size(timesteps));
%%

%Calculate A at time t using kinetics from previous timestep
%DL: check your prev/next timesteps...
for tt = 2:length(timesteps)
    
    dAdt = -kd(tt-1).*A(tt-1) + kf(tt-1).*(1-A(tt-1));
    dkfdt = (-kf(tt-1) + kf_0 + k_CamK.*Sigmoid(Ca(tt-1),(Ca_Kalpha.*(1-beta(tt-1))+Ca_beta.*beta(tt-1)),s_CamK))./tau_CamK;
    dkddt = (-kd(tt-1) + kd_0 + k_CaN.*Sigmoid(Ca(tt-1),Ca_CaN,s_CaN))./tau_CaN;
    dbetadt = (-beta(tt-1) + Sigmoid(Ca(tt-1),Ca_beta,s_beta))./tau_beta;
    
    A(tt) = A(tt-1) + dAdt.*dt;
    kf(tt) = kf(tt-1) + dkfdt.*dt;
    kd(tt) = kd(tt-1) + dkddt.*dt;
    beta(tt) = beta(tt-1) + dbetadt.*dt;
    
    Ca(tt) = Ca_0 + R(tt).*Ca_PSP.*A(tt);
end


%%
%viewdur = maxT;
timwin = [t_app-1000 maxT];
figure
subplot(5,1,1)
plot(timesteps,R,'k')
xlim(timwin)
ylabel('Spike Rate')
subplot(5,1,2)
plot(timesteps,A,'k')
ylabel('pGluA1')
xlim(timwin)
subplot(5,1,3)
plot(timesteps,Ca,'k')
ylabel('Ca')
xlim(timwin)
subplot(5,1,4)
plot(timesteps,kd,'b')
hold on
plot(timesteps,kf,'r')
ylabel('de/phos kinetics')
legend('kd','kf')
xlim(timwin)
subplot(5,1,5)
plot(timesteps,beta,'k')
ylabel('% Beta')
xlim(timwin)
xlabel('T (min)')

NiceSave('HomeostaticDynamics',figfolder,[],'includeDate',true)

