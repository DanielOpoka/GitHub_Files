function [] = WVD_plot(time, freq, fs,  ind_win, signal, FFT_norm, ...
    twin, fwin, tlim, flim, n_fontsize, title_label)

% % Time vs. Frequency
ax(3) = subplot(3, 3, [4, 8]);
hold on;
wvd(signal, fs, 'smoothedPseudo', twin, fwin)
hold off;
colormap jet;
colorbar off;

% Determining the order of magnitude for the time range
% tlim_ord = (-1 * floor(log10(tlim(2))));
% tmag_ord = 10^((floor(tlim_ord / 3)) * 3);
if time(end) < 1
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

    title_txt = "Smoothed Pseudo WVD of " + title_label;
    sgtitle(title_txt, 'FontSize', 30)
else
    xlim([tlim(1), tlim(2)]);
    ylim([flim(1), flim(2)] .* 1e-3);
    
    set(gca,'fontsize',n_fontsize);
%     Labelling
    ax(3).YTickLabel = ax(3).YTick .* 1e3;
    xlabel('Time [s]');
    ylabel('Frequency [Hz] ');
    grid on;

    title_txt = "Smoothed Pseudo WVD of " + title_label;
    sgtitle(title_txt, 'FontSize', 30)
end


signal = normalize(signal, 'norm', Inf);

% % Signal vs. Time
set(gcf, 'Position', get(0, 'Screensize'));
ax(1) = subplot(3, 3, [1, 2]);
hold on;
grid on;
% plot(time(ind_win(1):ind_win(2)) .* 1e9, twin, 'r')
plot(time(ind_win(1):ind_win(2)), twin, 'r')
% plot(time .* 1e9, signal, 'b')
plot(time, signal, 'b')
hold off;
% xlim([tlim(1), tlim(2)] .* 1e9);
xlim([tlim(1), tlim(2)]);
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
ax(2).YTickLabel = ax(3).YTick .* 1e3;
xlabel('FFT Amplitude [V/m]');

end



