

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200422';

%%

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);

%%
Run_SynHomeo(manip,parms,timeparms)

%% TTX Step: Increased Ca_PSP (TS)
timeparms.maxT = 72*60;
timeparms.preT = 5000;
parms.Ca_PSP = 2;
R_pre = 100;
R_post = 20;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);

%%
Run_SynHomeo(manip,parms,timeparms)

%% Excitability Pulse

R_baseline = 100;
R_pulse = [100 500];
t_pulse = [0 60];
pulsedur = 2;

manip.rate = @(t) R_baseline + ...
    R_pulse(1).*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse(2).*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur));

timeparms.maxT = 8*60;

Run_SynHomeo(manip,parms,timeparms)


%% Before/After Excitability 
parms= struct();
timeparms = struct();
timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 20;
t_TTX = 600;

%R_baseline = 100;
R_pulse = 800;
t_pulse = [-600 120 1800]+t_TTX;
pulsedur = 5;

manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX) + ...
    R_pulse.*(t>=t_pulse(1) & t<(t_pulse(1)+pulsedur)) + ...
    R_pulse.*(t>=t_pulse(2) & t<(t_pulse(2)+pulsedur)) + ...
    R_pulse.*(t>=t_pulse(3) & t<(t_pulse(3)+pulsedur));

%timeparms.maxT = 8*60;

Run_SynHomeo(manip,parms,timeparms)



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
