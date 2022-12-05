function [tau_RT] = ReverberationTime(time, r_t)

% Reverberation Time
num_tau_RT = []; dem_tau_RT = [];
for n_r_t = 1:size(r_t, 2)
    % Numerical Integration with trapz
    num_tau_RT_func = time.^2 .* r_t(:, n_r_t).^2;
    num_tau_RT(n_r_t) = trapz(time, num_tau_RT_func);

    dem_tau_RT_func = r_t(:, n_r_t).^2;
    dem_tau_RT(n_r_t) = trapz(time, dem_tau_RT_func);
end

tau_RT = sqrt(mean(num_tau_RT) / mean(dem_tau_RT)); % Reverberation Time [s]

end

