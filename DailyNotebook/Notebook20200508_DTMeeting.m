

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200508';

%% basic rate step (TTX)

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_WT',...
    'saveFig_AF',figfolder)

%% TTX Step: Increased Ca_PSP (TS)

parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;
parms.Ca_PSP = 2;
R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')

%% Block N: WT
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_WT')

%% Block M with TTX: WT
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockM = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockM_TTX_WT')
%% Block N: TS
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

parms.Ca_PSP = 2;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_TS');

%% Block N,M


parms= struct();
timeparms = struct();
manip = struct();
timeparms.maxT = 72*60;
timeparms.preT = 5000;

%    
%parms.kd_0 = 1e-4; 
%parms.kf_0 = parms.kd_0.*k0ratio;  

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);
manip.blockM = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

[simresults,simparms] = Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockNM_WT')

%%

Ca_Eq = -6.8;
b_eq = Sigmoid(Ca_Eq,simparms.Ca_beta,simparms.s_beta);
m_eq = Sigmoid(Ca_Eq,simparms.Ca_Kalpha-b_eq.*simparms.Ca_Kdelta,simparms.s_CamK);
n_eq = Sigmoid(Ca_Eq,simparms.Ca_CaN,simparms.s_CaN);
k0ratio = (m_eq.*simparms.k_CamK)./(n_eq.*simparms.k_CaN);

%% Excitability Pulse

R_baseline = 100;
R_pulse = [100 800];
t_pulse = [0 60];
pulsedur = 5;

manip.rate = @(t) R_baseline + ...
    R_pulse(1).*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse(2).*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur));

timeparms.maxT = 12*60;

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','Excitability_WT')


%% Before/After Excitability 
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 1200;

%R_baseline = 100;
R_pulse = 800;
t_pulse = [-1200 2400]+t_TTX;
pulsedur = 5;

manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX) + ...
    R_pulse.*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse.*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur));% + ...
  %  R_pulse.*(t>=t_pulse(3) & t<(t_pulse(3)+pulsedur));

%timeparms.maxT = 8*60;

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTX_Excitability')


%% Before/After Excitability : Beta Blocks
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 1200;

%R_baseline = 100;
R_pulse = 800;
t_pulse = [-1200 2400]+t_TTX;
pulsedur = 5;

manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX) + ...
    R_pulse.*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse.*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur));% + ...
  %  R_pulse.*(t>=t_pulse(3) & t<(t_pulse(3)+pulsedur));
  
manip.blockB = @(t) 0.*(t);

%timeparms.maxT = 8*60;

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTX_Excitability_BetaBlock')


%% basic rate step (TTX)

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 200;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BICStep_WT')



%% INternal/external calcium


R_baseline = 100;
R_pulse = [100 0 100];
Ca_pulse = [0 0.5 0.5];

t_pulse = [0];
pulsedur = 5;

for pp = 1:length(Ca_pulse)

manip.rate = @(t) R_baseline + ...
    R_pulse(pp).*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur));
manip.Ca_ext = @(t) 0 + ...
    Ca_pulse(pp).*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur));

timeparms.maxT = 12*60;

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname',['PulseExp',num2str(pp)])

end

%%

parms= struct();
timeparms = struct();
timeparms.maxT = 48*60;
timeparms.preT = 4000;

R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);

%%
Ca_PSPs = logspace(-0.5,0.7,48);
clear simresults_capsp
parfor cc = 1:length(Ca_PSPs)
    cc
    %bz_Counter(cc,length(Ca_PSPs),'sim')
    inparms = parms;
    inparms.Ca_PSP = Ca_PSPs(cc);
    simresults_capsp(cc) = Run_SynHomeo(manip,inparms,timeparms);
end

%%
AllSimResults = bz_CollapseStruct(simresults_capsp,1);

%%
figure
imagesc(AllSimResults.t_hr(1,:),log10(Ca_PSPs),AllSimResults.A)
hold on
plot([0 0],ylim(gca),'k')
LogScale('y',10)
xlim([-6 48])
ylabel('Ca_P_S_P');xlabel('t (hr)')
ColorbarWithAxis([0 0.25],'A')
axis xy

NiceSave('Vary_CaPSP',figfolder,[],'includeDate',true)
%%

R_posts = linspace(5,60,48);
clear simresults_rate
clear manip
parfor cc = 1:length(Ca_PSPs)
    cc
    %bz_Counter(cc,length(Ca_PSPs),'sim')
    inparms = parms;
    manip_varrate(cc).rate = @(t) R_posts(cc).*(t>=t_TTX) + R_pre.*(t<t_TTX);

    simresults_rate(cc) = Run_SynHomeo(manip_varrate(cc),inparms,timeparms);
end

%%
AllSimResults = bz_CollapseStruct(simresults_rate,1);

%%
figure
imagesc(AllSimResults.t_hr(1,:),(R_posts),AllSimResults.A)
hold on
plot([0 0],ylim(gca),'k')
%LogScale('y',10)
xlim([-6 48])
ylabel('PSP Rate: TTX');xlabel('t (hr)')
ColorbarWithAxis([0 0.25],'A')
NiceSave('Vary_RPost',figfolder,[],'includeDate',true)
