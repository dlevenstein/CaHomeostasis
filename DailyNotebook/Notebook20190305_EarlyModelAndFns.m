figfolder  = '/Users/dlevenstein/Project Repos/CaHomeostasis/DailyNotebook/Figures';
%%
%dA = -kf.*A + kd(1-A);
A = kf./(kf+fd);
Ca = Ca0 + (R_mini+R_spk).*Ca_PSP.*A;

kf = kf_0 + k_CamK.*Sigmoid(Ca,(Ca_Kalpha.*(1-beta)+Ca_beta.*beta),s_CamK);
kd = kd_0 + k_CaN.*Sigmoid(Ca,Ca_CaN,s_CaN);
beta = Sigmoid(Ca,Ca_beta,s_beta);


%% Parameters
kf_0 = 1;       %CamKII-independent (i.e. baseline) phosphorylation rate
kd_0 = 0.5;       %CaN-independent (i.e. baseline) dephosphorylation rate
k_CamK = 1;     %Maximal CamKII-mediated phoshorpylation rate
k_CaN = 2;      %Maximal CaN-mediated dephosphorylation rate

Ca_0 = 0;       %GluA1-independent (i.e. baseline) Calcium concentration    
Ca_PSP = 1;     %Calcium from 1Hz PSPs with GluA1 fully phosphorylated
Ca_Kalpha = 3.5;   %Ca midpoint for CamKII alpha
Ca_Kbeta = 2.5;  %Ca midpoint for CamKII beta
Ca_CaN = 2.25;     %Ca midpoint for CaN
Ca_beta = 0.5;    %Ca midpoint for CamKII alpha-beta transition

s_CamK = 10;    %Steepness of CamKII activation
s_CaN = 1;      %Steepness of CaN activation
s_beta = -4;     %Steepness of CamKII alpha-beta transition

tau_CamK = 1;   %Timescale of CamKII activation
tau_CaN = 10;   %Timescale of CaN activation
tau_beta = 100; %Timescale of CamKII alpha-beta transition


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

NiceSave('ActivationFunctions',figfolder,[],'includeDate',true)
