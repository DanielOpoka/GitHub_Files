function [W_alpha] = WignerDist(alpha, time, freq_WD)
% Wigner Distribution Function with alpha
W_alpha = zeros(length(time), length(freq_WD));

for ind_t = 1:length(time)
    for ind_f = 1:length(freq_WD)
        W_alpha(ind_t, ind_f) = (exp(2 * alpha * time(ind_t)) * ...
            (sin(2 * (2 * pi * freq_WD(ind_f)) * time(ind_t)) / (pi * (2 * pi * freq_WD(ind_f))))); 
    end
end
end
