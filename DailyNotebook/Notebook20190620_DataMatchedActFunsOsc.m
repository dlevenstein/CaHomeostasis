figfolder  = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Figures';


%% Parameters
kf_0 = 1;       %CamKII-independent (i.e. baseline) phosphorylation rate
kd_0 = 0.5;       %CaN-independent (i.e. baseline) dephosphorylation rate
k_CamK = 20;     %Maximal CamKII-mediated phoshorpylation rate
k_CaN = 20;      %Maximal CaN-mediated dephosphorylation rate

Ca_0 = -7.2;       %GluA1-independent (i.e. baseline) Calcium concentration
                    %(-6.854 at eqfrom SS)
Ca_PSP = 0.1;     %Calcium from 1Hz PSPs with GluA1 fully phosphorylated

Ca_Kalpha = -5.5;   %Ca midpoint for CamKII alpha (data)
%Ca_Kbeta = 0.5;  %Ca midpoint for CamKII beta
Ca_Kdelta = 1;  %How much CamKIIbeta lowers the "threshold" of CamK
                %10x more sensitive (data)
Ca_CaN = -6.2;     %Ca midpoint for CaN (data)
Ca_beta = -6.5;    %Ca midpoint for CamKII alpha-beta transition

s_CamK = 8;    %Steepness of CamKII activation (data)
s_CaN = 4;      %Steepness of CaN activation (data)
s_beta = -20;     %Steepness of CamKII alpha-beta transition

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
Ca_X = linspace(-7,-5,100);
R_X = linspace(0,20,100);
A_X = linspace(0,1,100);
[A_XY,R_XY] = meshgrid(A_X,R_X);


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

    
subplot(2,2,3)
    imagesc(A_X,R_X,Ca_0 + R_XY.*Ca_PSP.*A_XY)
    ColorbarWithAxis([-7 -5],'logCa')
    axis xy
    xlabel('A (pGluA1)');ylabel('Rate');
    title('Calcium Entry')
NiceSave('ActivationFunctions',figfolder,[],'includeDate',true)

%%
maxT = 7500;
dt = 0.01;
timesteps = 0:dt:maxT;

%%
R_pre = 40;
R_ttx = 10;
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
timwin = [t_app-2000 maxT];
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

