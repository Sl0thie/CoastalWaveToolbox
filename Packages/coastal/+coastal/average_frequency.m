function f = average_frequency(T)
%average_frequency Finds the average frequency of a wave
%   
    if(isnan(T));error('T is required');end
    f = 1 / T;

end

