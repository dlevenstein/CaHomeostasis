
%% Baseline Calcium
ca0s = {-7.75,-7.5,-7.25};
figure
for cc = 1:length(ca0s)
    Ca_0 = ca0s{cc};
%for Ca_psp = [0.01 0.1 1 10]

%Ca_0 = -7.5;
Ca_psp = 1;
kf_0 = 0.00001;
kd_0 = 0.00001;
kd_max = 0.3;
kf_max = 0.2;

m = 0;

s_n = 6;
s_b = -10;
s_m = 8;

Ca_m = -5.5;
Ca_n = -6.4;
Ca_b = -7;

Ca = linspace(Ca_0,-6.5,100);

%R = ((Ca - Ca_0).*(kf_0 + kd_0 + Sigmoid( Ca,Ca_n,s_n,kd_max)))./(kf_0.*Ca_psp);

n = Sigmoid( Ca,Ca_n,s_n,1);
kf = kf_0 + kf_max.*m;
kd = kd_0 + kd_max.*n;

A = kf./(kf+kd);



R = (Ca - Ca_0)./(A.*Ca_psp);


subplot(2,2,3)
hold on
plot(Ca,R)
xlabel('Ca');ylabel('R')
box off
ylim([0 200])
legend(cellfun(@num2str,ca0s,'UniformOutput',false),'location','northwest')
title('Ca_0')

subplot(2,2,1)
plot(Ca,Sigmoid( Ca,Ca_n,s_n,1))
hold on
plot(Ca,Sigmoid( Ca,Ca_b,s_b,1))
plot(Ca,Sigmoid( Ca,Ca_m,s_m,1))

subplot(2,2,4)
hold on
plot(Ca,A,'k')

end

%%

%% CaPSP
capsps = {0.01,0.1,0.5,1};
figure
for cc = 1:length(capsps)
    Ca_psp = capsps{cc};
%for Ca_psp = [0.01 0.1 1 10]

Ca_0 = -7.5;
%Ca_psp = 0.1;


Ca = linspace(Ca_0,-6,100);

%R = ((Ca - Ca_0).*(kf_0 + kd_0 + Sigmoid( Ca,Ca_n,s_n,kd_max)))./(kf_0.*Ca_psp);

n = Sigmoid( Ca,Ca_n,s_n,1);
kf = kf_0 + kf_max.*m;
kd = kd_0 + kd_max.*n;

A = kf./(kf+kd);



R = (Ca - Ca_0)./(A.*Ca_psp);


subplot(2,2,3)
hold on
plot(Ca,R)
xlabel('Ca');ylabel('R')
box off
ylim([0 200])
legend(cellfun(@num2str,capsps,'UniformOutput',false),'location','northwest')
title('Ca_p_s_p')

subplot(2,2,1)
plot(Ca,Sigmoid( Ca,Ca_n,s_n,1))
hold on
plot(Ca,Sigmoid( Ca,Ca_b,s_b,1))
plot(Ca,Sigmoid( Ca,Ca_m,s_m,1))

subplot(2,2,4)
hold on
plot(Ca,A,'k')

end



%% Log Ca! linear R->Ca relationship

capsps = {0.1,1,10,100};
figure
for cc = 1:length(capsps)
    Ca_psp = capsps{cc};
%for Ca_psp = [0.01 0.1 1 10]

Ca_0 = -7.5;
%Ca_psp = 0.1;


Ca = linspace(Ca_0,-6,100);

%R = ((Ca - Ca_0).*(kf_0 + kd_0 + Sigmoid( Ca,Ca_n,s_n,kd_max)))./(kf_0.*Ca_psp);

n = Sigmoid( Ca,Ca_n,s_n,1);
kf = kf_0 + kf_max.*m;
kd = kd_0 + kd_max.*n;

A = kf./(kf+kd);
R = 10.^(Ca - Ca_0)./(A.*Ca_psp);


subplot(2,2,3)
hold on
plot(Ca,R)
xlabel('Ca');ylabel('R')
box off
ylim([0 1000])
legend(cellfun(@num2str,capsps,'UniformOutput',false),'location','northwest')
title('Ca_p_s_p')

subplot(2,2,1)
plot(Ca,Sigmoid( Ca,Ca_n,s_n,1))
hold on
plot(Ca,Sigmoid( Ca,Ca_b,s_b,1))
plot(Ca,Sigmoid( Ca,Ca_m,s_m,1))

subplot(2,2,4)
hold on
plot(Ca,A,'k')

end

