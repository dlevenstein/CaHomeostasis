

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200706';

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
parms.Ca_PSP = 2;
Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')
%parms.Ca_PSP = 1; %(Back to default)


%% Block N: WT
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 24*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

%parms.kf_0 = 1e-4;

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_WT',...
    'FixYRange',true)

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

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockM_TTX_WT',...
    'FixYRange',true)
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

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockN_TS',...
    'FixYRange',true);

%% Block N,M


parms= struct();
timeparms = struct();
manip = struct();
timeparms.maxT = 24*60;
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

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','BlockNM_WT',...
    'FixYRange',true)