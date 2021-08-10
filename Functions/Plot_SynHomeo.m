function [] = Plot_SynHomeo(simresults,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'saveFig',false)
addParameter(p,'figname',[])
addParameter(p,'FixYRange',false)
addParameter(p,'YScale',1)
addParameter(p,'figwidth',1)
addParameter(p,'fignum',1)
addParameter(p,'title',[])
addParameter(p,'linecolor','k')


parse(p,varargin{:})
saveFig = p.Results.saveFig;
figname = p.Results.figname;
FixYRange = p.Results.FixYRange;
figwidth = p.Results.figwidth;
fignum = p.Results.fignum;
plottitle = p.Results.title;
linecolor = p.Results.linecolor;

%%
timwin = [-5 simresults.t_hr(end)];

subplot(6,figwidth,(0.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.R,'color',linecolor,'linewidth',2)
    if isfield(simresults,'p')
       % plot(simresults.t_hr,
    end
    xlim(timwin)
    box off
    
    ylabel('Spike Rate')
    if FixYRange
        ylim([0 150])
    end
    title(plottitle)
subplot(6,figwidth,(1.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.A,'color',linecolor,'linewidth',2)
    ylabel('pGluA1')
    xlim(timwin)
    %ylim([0 1])
    box off
	if FixYRange
        ylim([0 0.2].*FixYRange)
    end
subplot(6,figwidth,(2.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.Ca,'color',linecolor,'linewidth',2)
    ylabel('Ca')
    xlim(timwin)
    box off
    if FixYRange
        ylim([-8 -6])
    end
    %ylim([Ca_0 -2.5])
subplot(6,figwidth,(3.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.m,'r','linewidth',2)
    hold on
    ylabel('CamK Gate')
    xlim(timwin)
    box off
    if FixYRange
        ylim([0 0.015].*FixYRange)
    end
    %ylim([0 1])
    
subplot(6,figwidth,(4.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.n,'b','linewidth',2)
    ylabel('CaN Gate')
    xlim(timwin)
    box off
    if FixYRange
        ylim([0 0.1].*FixYRange)
    end
subplot(6,figwidth,(5.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.b,'color',linecolor,'linewidth',2)
    ylabel('% Beta')
    xlim(timwin)
    box off
    xlabel('T (hr)')
    %ylim([0 1])
    if FixYRange
        ylim([0 1])
    end

if saveFig  
    NiceSave(figname,saveFig,[],'includeDate',true)
end

end

