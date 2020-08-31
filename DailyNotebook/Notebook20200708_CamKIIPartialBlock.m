

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200708';



%% Block M with TTX: WT
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

blocks = linspace(0,1,11);


for bb = 1:length(blocks)
    bz_Counter(bb,length(blocks),'Sim')
        manip.blockM = @(t) blocks(bb).*(t>=t_TTX) + 1.*(t<t_TTX);

        simresults(bb) = Run_SynHomeo(manip,parms,timeparms);
end
%%
simresults_combined = bz_CollapseStruct(simresults,1);
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

Plot_SynHomeo(simresults(11),'FixYRange',1.75,'figwidth',2,'fignum',2)
Plot_SynHomeo(simresults(2),'FixYRange',1.75,'figwidth',2,'fignum',2)
Plot_SynHomeo(simresults(1),'FixYRange',1.75,'figwidth',2,'fignum',2,'title','0, 90, 100%')

NiceSave('PartialBlockN',figfolder,[],'figtype','pdf') 

%%

