function [Ca] = CaFromA(A,R,Ca_0,Ca_PSP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%A = (kf_0 + k_CamK.*m)./(kf_0 + k_CamK.*m + kd_0 + k_CaN.*n);

Ca = Ca_0 + R.*Ca_PSP.*A;

end

