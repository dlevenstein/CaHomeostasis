

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200731';

%% basic rate step (TTX)

%manip.rate = @(t) 20.*ones(size(t));

parms= struct();
timeparms = struct();
manip = struct();

parms.kf_0 = 1e-4;      %(default: 1e-4) %CamKII-independent (i.e. baseline) phosphorylation rate
parms.kd_0 = 1e-7;  %(default: 1e-5)
parms.k_CaN = 0.11; %(default: 0.15)
parms.Ca_beta = -7.0; %(default: -7.05)
parms.s_beta = -35; %(default: -35)
parms.k_CamK = 0.15; %default: 0.1

timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 150;
R_post = 25;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_WT')


% TTX Step: Increased Ca_PSP (TS)
% parms.Ca_PSP = 2;
% Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')
%parms.Ca_PSP = 1; %(Back to default)
%% Early phase: tau_b and kf_o
timeparms.maxT = 36*60;
taus_b = [100 300 500];

clear simresults
for bb = 1:length(taus_b)
    parms.tau_beta = taus_b(bb);
    simresults_taubeta(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
end
 parms.tau_beta = 300; %Back to default
 
 kfs_0 = [0.5e-4 1e-4 2e-4];

for bb = 1:length(kfs_0)
    parms.kf_0 = kfs_0(bb);
    simresults_kf0(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
end
parms.kf_0 = 1e-4;
%%
%results_taubeta = bz_CollapseStruct(simresults,1,'justcat',true);

%%
figure
for bb = 1:length(taus_b)
    Plot_SynHomeo(simresults_taubeta(bb),'FixYRange',1,'figwidth',2,'title','tau_b')
end

for bb = 1:length(taus_b)
    Plot_SynHomeo(simresults_kf0(bb),'FixYRange',1.75,'figwidth',2,'fignum',2,'title','kf_0')
end
NiceSave('EarlyPhase',figfolder,[]) 
%% Peak Magnitude: k_cam
kms = [0.15 0.225 0.3 0.4];

clear simresults
for bb = 1:length(kms)
    parms.k_CamK = kms(bb);
    simresults_km(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
end
 parms.k_CamK = 0.1; %Back to default
 
%  kfs_0 = [0.5e-4 1e-4 2e-4];
% 
% for bb = 1:length(kfs_0)
%     parms.kf_0 = kfs_0(bb);
%     simresults_kf0(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
% end
% parms.kf_0 = 1e-4;

%%
figure
for bb = 1:3
    Plot_SynHomeo(simresults_km(bb),'FixYRange',false,'figwidth',2,'title','k_C_a_m_K')
end
Plot_SynHomeo(simresults_km(4),'FixYRange',false,'figwidth',2,'fignum',2,'title','High k_C_a_m_K')

NiceSave('CamKPeak',figfolder,[]) 
%% Peak Magnitude: k_caN
kns = [0.05 0.1 0.15];

clear simresults
for bb = 1:length(kns)
    parms.k_CaN = kns(bb);
    simresults_kn(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
end
 parms.k_CaN = 0.11; %Back to default
 
%  kfs_0 = [0.5e-4 1e-4 2e-4];
% 
% for bb = 1:length(kfs_0)
%     parms.kf_0 = kfs_0(bb);
%     simresults_kf0(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
% end
% parms.kf_0 = 1e-4;

%%
figure
for bb = 1:3
    Plot_SynHomeo(simresults_kn(bb),'FixYRange',false,'figwidth',2,'title','k_C_a_N')
end

NiceSave('DampedOsc',figfolder,[]) 


%% Tau_b and Kf_0 for bigger peak at same time

parms.tau_beta = taus_b(1);
parms.kf_0 = kfs_0(1);
simresults_net(bb) = Run_SynHomeo(manip,parms,timeparms,'showfig',false);
 parms.tau_beta = 300; %Back to default
 parms.kf_0 = 1e-4;
%%
figure
Plot_SynHomeo(simresults_net(bb),'FixYRange',false,'figwidth',2,'title','both')