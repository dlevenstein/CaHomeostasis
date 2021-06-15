

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210119';



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


[simresults,simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','TTXStep_WT','showfig_ActFun',true,...
    'saveFig_AF',figfolder)

%%
%Earlier peak?
taus_pre = [300 600 1200];
for pp = 1:length(taus_pre)
parms.tau_pre = taus_pre(pp);
[simresults(pp),simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname',['SlowTauPre_',num2str(taus_pre(pp))])
parms.tau_beta = taus_pre(pp);
[simresults(pp),simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname',['SlowTauNuc_',num2str(taus_pre(pp))])
end

%%
parms.tau_pre = 600;
parms.tau_beta = 600;
parms.kf_0 = 1e-4;
parms.k_CamK = 0.15;

parms.kf_0 = 2e-4;
parms.k_CamK = 0.30;
[simresults(pp),simparms] = Run_SynHomeo_Pre(manip,parms,timeparms,...
    'saveFig',figfolder,'figname','PhosPeak')
%%
figure
% subplot(2,2,1)
% plot(simresults_combined.t_hr(1,:),simresults_combined.A([1:end],:))
% xlim([-6 24])
% ylim([0 0.2])
% xlabel('t (hr)');ylabel('A')
% for bb = 1:length(blocks)
%     Plot_SynHomeo(simresults(bb),'FixYRange',1.75,'figwidth',2,'fignum',1,'title','0-100%')
% end

Plot_SynHomeo(simresults,'FixYRange',1.75,'figwidth',2,'fignum',2)
%NiceSave('PartialBlockN',figfolder,[],'figtype','pdf') 

%%

