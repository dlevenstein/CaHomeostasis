figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200701';
%%
Ca_0 = -7.5;       %GluA1-independent (i.e. baseline) Calcium concentration
                            %(-6.854 at eqfrom SS)
Ca_PSP = 1;     %Calcium per 1Hz PSPs with GluA1 fully phosphorylated
kf_0 = 1e-4;       %CamKII-independent (i.e. baseline) phosphorylation rate
Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)

Ca_B = -7; %Baseline Ca
R_B = logspace(0,3.5,22); %Baseline Rate

A_B = (10.^(Ca_B-Ca_0)./(R_B.*Ca_PSP));
n_B = Sigmoid(Ca_B,Ca_CaN,s_CaN);

%CaN dephos rate to give baseline Ca at baseline Rate
kN = kf_0.*(1-A_B)./(A_B.*n_B)


%% basic rate step (TTX)

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();



timeparms.maxT = 0;
timeparms.preT = 4000;

for nn = 1:length(kN)
    bz_Counter(nn,length(kN),'sim')
parms.k_CaN = kN(nn);

R_pre = R_B(nn);
manip.rate = @(t) R_pre .* ones(size(t));

simresults = Run_SynHomeo(manip,parms,timeparms);
equilCa_2(nn) = simresults.Ca(end);
close all
end
%%

kNforRate = kN([10 13 17])
R_pres = logspace(0,3.5,22);

for nn = 1:length(kNforRate)
parms.k_CaN = kNforRate(nn);
for rr = 1:length(R_pres)
    bz_Counter(rr,length(R_pres),'sim')
    
R_pre = R_pres(rr);
manip.rate = @(t) R_pre .* ones(size(t));

simresults = Run_SynHomeo(manip,parms,timeparms);
equilCa(rr,nn) = simresults.Ca(end);
close all
end

end


%%
figure
subplot(2,2,2)
plot(log10(R_pres),equilCa)
hold on
plot(xlim(gca),[1 1].*Ca_B,'k--')
LogScale('x',10)
xlabel('Rate')
ylabel('Equilibrium Ca')
legend(num2str(kN(10)), num2str(kN(13)), num2str(kN(17)),'location','northwest')

subplot(2,2,1)
plot(log10(R_B(kN>=0)),log10(kN(kN>=0)),'k')
hold on
plot(log10(R_B([10 13 17])),log10(kN([10 13 17])),'o')
xlabel('Baseline Rate');ylabel(['k_C_a_N for Ca = ',num2str(Ca_B)])
LogScale('xy',10)

% subplot(2,2,3)
% plot(log10(R_B),equilCa_2)
% hold on
% plot(xlim(gca),[1 1].*Ca_B,'k--')
% LogScale('x',10)
% xlabel('Baseline Rate');ylabel('Ca Eq')

NiceSave('BaselineCa',figfolder,[],'includeDate',true)
