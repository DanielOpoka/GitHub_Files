function [] = STFT_plot(time, freq, fs,  ind_win, win, signal, FFT_norm, ...
    over_len, tlim, flim, n_fontsize, title_label)

% % Signal vs. Time
set(gcf, 'Position', get(0, 'Screensize'));
ax(1) = subplot(3, 3, [1, 2]);
hold on;
grid on;
plot(time(ind_win(1):ind_win(2)) .* 1e9, win, 'r')
plot(time .* 1e9, signal, 'b')
hold off;
xlim([tlim(1), tlim(2)] .* 1e9);
set(gca,'fontsize',n_fontsize);
% Labelling
ylabel('Normalized Amplitude [V/m]');

% % FFT vs. Freq
ax(2) = subplot(3, 3, [6, 9]);
hold on;
grid on;
plot(FFT_norm, freq .* 1e-9,'b');
hold off;
ylim([flim(1), flim(2)] .* 1e-9);
set(gca,'fontsize',n_fontsize);
% Labelling
xlabel('FFT Amplitude [V/m]');

% % Time vs. Frequency
ax(3) = subplot(3, 3, [4, 8]);
hold on;
stft(signal, fs, 'Window', win, 'OverlapLength', over_len, ... 
    'FrequencyRange', 'onesided', 'FFTLength', length(time));
hold off;
colormap jet;
colorbar off;

% Determining the order of magnitude for the time range
% tlim_ord = (-1 * floor(log10(tlim(2))));
% tmag_ord = 10^((floor(tlim_ord / 3)) * 3);
if (time(end) > 1e-6)
    tmag_ord = 1e6;
    xlim([tlim(1), tlim(2)] .* tmag_ord);
    ylim([flim(1), flim(2)] .* 1e-9);
    ax(3).XTickLabel = ax(3).XTick .* 1000;
else
    tmag_ord = 1e9;
    xlim([tlim(1), tlim(2)] .* tmag_ord);
    ylim([flim(1), flim(2)] .* 1e-9);
end

set(gca,'fontsize',n_fontsize);
% Labelling
% ax(3).XTickLabel = ax(3).XTick .* 1000;
xlabel('Time [ns]');
ylabel('Frequency [GHz] ');
grid on;

title_txt = "STFT of " + title_label;
sgtitle(title_txt, 'FontSize', 30)

end

