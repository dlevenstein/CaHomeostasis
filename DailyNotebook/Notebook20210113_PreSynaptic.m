

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20210113';



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

%% Kinetic (Phosphorylation) equation
kd0 = 1e-2;
kdm = 2e-2;
kf0 = 1e-2;
kfm = 2e-2;

Ca = linspace(simparms.Ca_0,-5,100);
n = Sigmoid(Ca,simparms.Ca_CaN,simparms.s_CaN);

kd = kdm.*n + kd0;
kf = kfm.*(Ca-simparms.Ca_0) + kf0;

S_inf = kd./(kf+kd);
tau_s = 1./(kd+kf);

figure
subplot(2,2,1)
    plot(n,kd,'linewidth',2)
    xlabel('n');ylabel('Presynaptic kd')
subplot(2,2,2)
    plot(Ca,kf,'linewidth',2)
    hold on
    plot(Ca,kd,'linewidth',2)
    xlabel('Ca');ylabel('Presynaptic k')
    legend('kf','kd','location','northwest')
subplot(2,2,3)
    plot(Ca,S_inf,'linewidth',2)
    xlabel('Ca');ylabel('S inf (steady state)')
subplot(2,2,4)
    plot(Ca,tau_s,'linewidth',2)
    xlabel('Ca');ylabel('tao_s')

NiceSave('S_curves',figfolder,'Presyn_Phos')
%%
kd = kdm.*simresults.n + kd0;
kf = kfm.*(simresults.Ca-simparms.Ca_0) + kf0;

dt = mode(diff(simresults.t_sec));
s = zeros(size(simresults.t_sec));
for tt = 2:length(simresults.t_sec)   
    dsdt = kd(tt-1).*(1-s(tt-1)) - kf(tt-1).*s(tt-1);
    s(tt) = s(tt-1) + dsdt.*dt;
end
%%
figure
subplot(4,1,1)
plot(simresults.t_hr,kd)
subplot(4,1,2)
plot(simresults.t_hr,kf)
subplot(4,1,4)
plot(simresults.t_hr,s)

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

