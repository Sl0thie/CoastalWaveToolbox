function OMEGA = angular_frequency(T)
%angular_frequency Summary of this function goes here
%   Detailed explanation goes here
    if(isnan(T));error('T is required');end
    OMEGA = (pi * 2) / T;

end

