

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210402';



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

R_pre = 100;
R_post = 0;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);


parms.Ca_pre = -7;
parms.mini_max = 80;

[simresults,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','TTXStep_WT_Bell','showfig_ActFun',true,...
    'saveFig_AF',figfolder,'PreActFun','bell')

%% TTX Step: Increased Ca_PSP (TS)
parms.Ca_PSP = 1.2;
parms.Ca_pre = -7;
parms.mini_max = 80;
parms.Ca_Kalpha = -5.5;
%parms.tau_pre = 500;
%parms.tau_beta = 500;
[simresults,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','TTXStep_TS_Bell','showfig_ActFun',true,...
    'saveFig_AF',figfolder,'PreActFun','bell');


%% Block N: WT


timeparms.maxT = 72*60;
timeparms.preT = 5000;

R_pre = 100;
R_post = 100;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.blockN = @(t) 0.*(t>=t_TTX) + 1.*(t<t_TTX);

[simresults,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','TTXStep_BlockN_WT_Bell','showfig_ActFun',true,...
    'saveFig_AF',figfolder,'PreActFun','bell')





