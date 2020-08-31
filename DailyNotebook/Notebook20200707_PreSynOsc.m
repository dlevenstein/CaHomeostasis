

figfolder  = '/Users/dl2820/Project Repos/CaHomeostasis/DailyNotebook/Notebook20200707';

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

%% Presynaptic frequency fluctuation!
R_res = @(t,R_eq,A,tau,f) R_eq-A.*exp(-t./tau).*cos((t).*2*pi./f);

R_eq = 30; %Equilibrium rate
A = 20;     %Immediate post-TTX rate
tau = 36.*60;   %Decay time constant
f = 36.*60;     %Presynaptic oscillation frequency
%phi = 2;


%manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip.rate = @(t) R_res(t,R_eq,A,tau,f).*(t>t_TTX) + R_pre.*(t<t_TTX);


time = linspace(-6,72,1000).*60;
figure
plot(time,manip.rate(time),'k')
box off
xlabel('Time (min)')
ylabel('PSP Rate')

%%

Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_WT')


% TTX Step: Increased Ca_PSP (TS)
parms.Ca_PSP = 2;
Run_SynHomeo(manip,parms,timeparms,'saveFig',figfolder,'figname','TTXStep_TS')
%parms.Ca_PSP = 1; %(Back to default)


