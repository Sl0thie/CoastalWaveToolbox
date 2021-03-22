function k = wave_number(L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if(isnan(L));error('L is required');end
    k = (pi * 2) / L;
    
end

