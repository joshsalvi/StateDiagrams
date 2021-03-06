function [G,Xc,ramp_time] = space_calc(F_min,F_max,N_F,k_min,k_max,N_k,Xc_slope_max,G_slope_max,alpha,beta,k_sf)
%Area_calc calculates the state space in terms of gain and command
%displacement
%
%[G,Xc,ramp_time] = space_calc(F_min,F_max,N_F,k_min,k_max,N_k,F_slope_max,G_slope_max,alpha,beta,k_sf)
%
%G: two-dimensional matrix with values for gain
%Xc: two-dimensional matrix with values for command displacement
%ramp_time: time for all ramps in the protocol
%
%F_min,F_max: values for minimum and maximum force, respectively
%k_min,k_max: values for minimum and maximum stiffness, respectively
%N_F,N_k: number of points in force and stiffness
%Xc_slope_max,G_slope_max: maximum slope for gain and offset ramps
%alpha,beta: calibration parameters
%k_sf: stimulus fiber stiffness
%
%If N_k or N_F == 1, then this script will assign F=F_min and/or k=k_min;
%
%All values are required.
if N_F > 1
    if F_min == F_max
        F(1:N_F) = F_min;
    else
        F = F_min:((F_max-F_min)/(N_F-1)):F_max;
    end
    
else
    F = F_min;
end

if N_k > 1
    if k_min == k_max
        k(1:N_k) = k_min;
    else
        k = k_min:((k_max-k_min)/(N_k-1)):k_max;
    end
else
    k = k_min;
end
% Calculate space
% Map points in terms of G and Xc
for i = 1:N_k
    for j = 1:N_F
        G(i,j) = (k(i)-k_sf)/(alpha*beta*k_sf);
        Xc(i,j) = F(j)/(alpha*beta*G(i,j)*k_sf);
    end
end

% Determine maximum ramp times for G and Xc
ramp_G = max(abs(G))/G_slope_max;
ramp_Xc = max(abs(Xc))/Xc_slope_max;

if ramp_G > ramp_Xc
    ramp_time = ramp_G;
else
    ramp_time = ramp_Xc;
end
