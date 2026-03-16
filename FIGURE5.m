load('v1');
load('a_resampled');

figure('Position',[100,100,1400,300])
time = (1:length(v1))/50 ;  
plot(time,a_resampled/max(abs(a_resampled))*1.5,'LineWidth',2.5,'Color',  [0, 0.4470, 0.7410]);
hold on

line1 = plot(time,v1/max(abs(v1)),'-.r','LineWidth',2.5);
ylim([-1.5 1.5]);
xlim([time(1),time(end)])
xlabel('Time (s)')
ylabel('Normalized waveform');
set(gca ,'Fontsize',13);
legend('Photonic NLOS','GT');

load('net_trained.mat');
load('feature.mat');
load('long_signal.mat');
fs=10;window_sec =150; window_pts = round(window_sec * fs); 

step = fs * 1; 
num_steps = floor((length(long_signal) - window_pts) / step);

time_axis = (0:num_steps-1) * step / fs;
probs_store = zeros(num_steps, 5); 
hWait = waitbar(0, 'Predicting...');

for i = 1:num_steps
    if mod(i, 200) == 0; waitbar(i/num_steps, hWait); end
    
    start_pt = (i-1)*step + 1;
    end_pt = start_pt + window_pts - 1;
    
    current_window = long_signal(start_pt:end_pt);
    
    % 
    sig_input = (current_window - mean(current_window)) / (std(current_window) + 1e-6);
    % Reshape [600, 1, 1, 1]
    dl_wave = dlarray(reshape(sig_input, [window_pts, 1, 1, 1]), 'SSCB');
    
    % 
    feats = extract_features(current_window, fs);
    % 
    feats_norm = (feats - feat_mean) ./ (feat_std + 1e-6);
    dl_feat = dlarray(feats_norm', 'CB');
    
    % 
    if canUseGPU
        dl_wave = gpuArray(dl_wave);
        dl_feat = gpuArray(dl_feat);
    end
    
    dl_out = predict(net, dl_wave, dl_feat); % 
    scores = extractdata(gather(dl_out));
    
    probs_store(i, :) = scores';
end
close(hWait);
% 6. 
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

ax1 = subplot(3, 1, 1);
t_full = (0:length(long_signal)-1)/fs;
plot(t_full, long_signal, 'Color', [0.2, 0.2, 0.2]); hold on;
yl = ylim;
% 
fill([3000, 4500, 4500, 3000]/10, [yl(1) yl(1) yl(2) yl(2)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([6000, 7500, 7500, 6000]/10, [yl(1) yl(1) yl(2) yl(2)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([9000, 10740, 10740,9000]/10, [yl(1) yl(1) yl(2) yl(2)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.1, 'EdgeColor', 'none');   
fill([12500, 14500, 14500, 12500]/10, [yl(1) yl(1) yl(2) yl(2)], [0.8500 0.3250 0.0980], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('Real time respiratory waveform detection', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized amplitude'); grid on; xlim([0, length(long_signal)/fs]);
xlim([75 1450-75]);

ax2 = subplot(3, 1, 3);
plot_time = time_axis+ window_sec / 2;
h=area(plot_time, probs_store);
title('Classification results', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Confidence level'); ylim([0 1]); grid on;
new_order_names = {'Normal', 'Fast', 'Slow', 'CS', 'Ataxic'};
% 3. 
new_order_idx = [1, 4, 5, 2, 3];
% 4. 
legend(h(new_order_idx), new_order_names, 'Location', 'EastOutside');
 xlim([75 1450-75]);xlabel('Time (s)');

ax3 = subplot(3, 1, 2);
is_abnormal = probs_store(:, 1) < 0.5; % 
is_abnormal_smooth = smoothdata(double(is_abnormal), 'movmean', 5);
plot(plot_time, is_abnormal_smooth, 'k', 'LineWidth', 2);
area(plot_time, is_abnormal_smooth, 'FaceColor', 'r', 'FaceAlpha', 0.3);
title('Alarm status', 'FontSize', 12, 'FontWeight', 'bold');
 yticks([0 1]); yticklabels({'Normal', 'Abnormal'});
ylim([-0.1 1.5]); grid on;
linkaxes([ax1, ax2, ax3], 'x'); xlim([0, length(long_signal)/fs]);
xlim([75 1450-75]);

function feats = extract_features(sig, fs)
    % 
    sig = sig(:); 
    
    % --- 1
    sig = detrend(sig);
    if std(sig) > 1e-6
        sig = (sig - mean(sig)) / std(sig); 
    else
        sig = zeros(size(sig)); 
    end
    
    % --- 2. 
    [Pxx, F] = periodogram(sig, [], 1024, fs);
    [~, max_idx] = max(Pxx); 
    f_dom_freq = F(max_idx); % 
    
    % 
    p_norm = Pxx / (sum(Pxx) + 1e-10);
    f_entropy = -sum(p_norm .* log2(p_norm + 1e-10));

    % --- 3. 
    [up_env, ~] = envelope(sig, round(fs*1.5), 'rms'); % 
    
    % 
    f_env_cv = std(up_env) / (mean(up_env) + 1e-6);
    
    % 
    [pks, ~] = findpeaks(sig, 'MinPeakDistance', fs*2); % 
    if length(pks) > 3
        f_peak_cv = std(pks) / (mean(pks) + 1e-6); % 
    else
        f_peak_cv = 1.0; %
    end

    % --- 4. 
    win_len = round(fs * 4); 
    sig_std_local = movstd(sig, win_len);
    adaptive_thresh = 0.35; % 
    is_silent = sig_std_local < adaptive_thresh;
    f_silence_ratio = sum(is_silent) / length(sig); % 

    % 
    diff_silent = diff([0; is_silent(:); 0]);
    starts = find(diff_silent == 1);
    ends = find(diff_silent == -1);
    durations = (ends - starts) / fs;
    if length(durations) >= 2
        f_silence_reg = f_silence_ratio*std(durations) / (mean(durations) + 1e-6);
    else
        f_silence_reg = (f_silence_ratio > 0.1) * 1.0; % 
    end

    % --- 5.
    if length(up_env) > fs*10
        [xc, lags] = xcorr(up_env - mean(up_env), 'coeff');
        half_xc = xc(lags > 0);
        f_env_corr = max(half_xc); % 
    else
        f_env_corr = 0;
    end

    
    feats = [f_dom_freq, f_entropy, f_env_cv, f_peak_cv, f_silence_ratio, f_silence_reg, f_env_corr];
end

%%

close all;
figure('Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 20, 20]); 


C_GT = [0.85, 0.16, 0.0];      % Ground Truth 
C_Main = [0, 0.447, 0.741];    % Sensor 
C_Wave = [0.15, 0.15, 0.15];   % 
Colors = [
    0.30, 0.75, 0.30; % Normal - 
    0.85, 0.15, 0.15; % Fast - 
    0.15, 0.45, 0.70; % Slow - 
    0.95, 0.60, 0.10; % CS - 
    0.45, 0.25, 0.60  % Ataxic - 
];


% ax1 = subplot(4, 1, 1);
% t_v1 = (1:length(v1))/50;
% p1 = plot(t_v1, a_resampled/max(abs(a_resampled))*1.1, 'Color', C_Main, 'LineWidth', 1.5); hold on;
% p2 = plot(t_v1, v1/max(abs(v1)), '-.', 'Color', C_GT, 'LineWidth', 1.2);
% ylabel({'Norm.','Wavef.'});
% xlabel('Time (s)'); % 
% title('a', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
% legend([p1, p2], {'Photonic NLOS', 'Reference (GT)'}, 'Box', 'off', 'Orientation', 'horizontal', 'Location', 'northoutside');
% grid on; xlim([0, 60]); ylim([-1.5, 1.5]);
% set(ax1, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial');
% 
% % 
% text(0.98, 0.1, sprintf('Pearson r = 0.784'), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
% 
% 
 t_crop = [75, 1375]; % 
ax1 = axes('Position', [0.13, 0.78, 0.775, 0.18]); % 
try
    img = imread('AA1.png'); % 
    imshow(img, 'Parent', ax1);
    axis(ax1, 'off'); % 
    text(ax1, -0.08, 1.1, 'a', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
catch
    warning('AA.png not found. Skipping panel a.');
    delete(ax1); % 
end

ax2 = subplot(4, 1, 2);
t_full = (0:length(long_signal)-1)/fs;
plot(t_full, -long_signal, 'Color', C_Wave, 'LineWidth', 0.8); hold on;
yl = ylim;
fill_ranges = {[300, 450], [600, 750], [900, 1074], [1250, 1450]};
for i = 1:length(fill_ranges)
    f = fill([fill_ranges{i}(1) fill_ranges{i}(2) fill_ranges{i}(2) fill_ranges{i}(1)], ...
             [yl(1) yl(1) yl(2) yl(2)], [0, 0, 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
ylabel({'Norm.','Wavef.'});
set(gca, 'XTickLabel', [], 'XColor', 'none');
title('b', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
grid on; axis tight; xlim(t_crop);


ax3 = subplot(4, 1, 3);
is_abnormal_smooth = smoothdata(double(probs_store(:, 1) < 0.5), 'movmean', 5);
area(plot_time, is_abnormal_smooth, 'FaceColor', [0.8, 0.2, 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.6, 0, 0]);
ylabel({'Alarm','Status'});
yticks([0 1]); yticklabels({'Normal', 'Abnormal'});
set(gca, 'XTickLabel', [], 'XColor', 'none');
title('c', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
grid on; xlim(t_crop); ylim([-0.1 1.5]);


ax4 = subplot(4, 1, 4);
h = area(plot_time, probs_store, 'LineStyle', 'none');
for k = 1:5, h(k).FaceColor = Colors(k,:); h(k).FaceAlpha = 0.85; end
ylabel({'Confid.','Level'});
xlabel('Time (s)');
title('d', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
lgd = legend(h(new_order_idx), new_order_names, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
grid on; xlim(t_crop); ylim([0 1]);

linkaxes([ax2, ax3, ax4], 'x'); %
set([ax2, ax3, ax4], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);

pos1 = get(ax1, 'Position');
set(ax1, 'Position', [pos1(1), pos1(2)+0.02, pos1(3), pos1(4)]);

%
close all;
close all;
figure('Color', 'w', 'Units', 'centimeters', 'Position', [1, 1, 25, 20]); 

left = 0.15;      % 
width = 0.78;     % 
h1 = 0.33;  b1 = 0.6% 
h2 = 0.12;  b2 = 0.43;  % 
h3 = 0.12;  b3 = 0.28;  % 
h4 = 0.08;  b4 = 0.15; % 
tX = -0.0;       % 
tY = 1.2;        % 
tFS = 12;         % 
tFW = 'bold';     % 
tFN = 'Arial';    % 
C_GT = [0.85, 0.16, 0.0]; C_Main = [0, 0.447, 0.741]; C_Wave = [0.15, 0.15, 0.15];
Colors = [0.3,0.75,0.3; 0.85,0.15,0.15; 0.15,0.45,0.7; 0.95,0.6,0.1; 0.45,0.25,0.6];

% ================= 
ax1 = axes('Position', [left, b1, width, h1]);
try
    img = imread('AA1.png'); 
    imshow(img, 'Parent', ax1);
    axis(ax1, 'off');
    text(ax1, tX, tY, 'a', 'Units', 'normalized', 'FontSize', tFS, 'FontWeight', tFW, 'FontName', tFN);
catch
    text(0.5, 0.5, 'Image Not Found', 'HorizontalAlignment', 'center');
end

% ================= 
ax2 = axes('Position', [left, b2, width, h2]);
t_full = (0:length(long_signal)-1)/fs;
plot(t_full, -long_signal, 'Color', C_Wave, 'LineWidth', 0.8); hold on;
yl = ylim;
fill_ranges = {[300, 450], [600, 750], [900, 1074], [1250, 1450]};
for i = 1:length(fill_ranges)
    fill([fill_ranges{i}(1) fill_ranges{i}(2) fill_ranges{i}(2) fill_ranges{i}(1)], ...
         [yl(1) yl(1) yl(2) yl(2)], [0, 0, 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
end
ylabel({'Norm.','Wavef.'});
set(ax2, 'XTickLabel', [], 'XColor', 'none');
text(ax2, tX, tY, 'b', 'Units', 'normalized', 'FontSize', tFS, 'FontWeight', tFW, 'FontName', tFN);
%title('b', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
grid on; xlim(t_crop);

% ================= 
ax3 = axes('Position', [left, b3, width, h3]);
is_abnormal_smooth = smoothdata(double(probs_store(:, 1) < 0.5), 'movmean', 5);
area(plot_time, is_abnormal_smooth, 'FaceColor', [0.8, 0.2, 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.6, 0, 0]);
ylabel({'Alarm','Status'});
yticks([0 1]); %yticklabels({'Normal', 'Abnormal'});
set(ax3, 'XTickLabel', [], 'XColor', 'none');
text(ax3, tX, tY, 'c', 'Units', 'normalized', 'FontSize', tFS, 'FontWeight', tFW, 'FontName', tFN);
grid on; xlim(t_crop); ylim([-0.1 1.5]);

ax4 = axes('Position', [left, b4, width, h4]);
h = area(plot_time, probs_store, 'LineStyle', 'none');
for k = 1:5, h(k).FaceColor = Colors(k,:); h(k).FaceAlpha = 0.85; end
ylabel({'Confid.','Level'});
xlabel('Time (s)');
text(ax4, tX, tY, 'd', 'Units', 'normalized', 'FontSize', tFS, 'FontWeight', tFW, 'FontName', tFN);
grid on; xlim(t_crop); ylim([0 1]);

linkaxes([ax2, ax3, ax4], 'x'); 
set([ax2, ax3, ax4], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);

lgd = legend(ax4, h(new_order_idx), new_order_names, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
lgd_pos = get(lgd, 'Position');
lgd_pos(2) = lgd_pos(2) - 0.05; 
set(lgd, 'Position', lgd_pos);
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLegendFontName', 'Times New Roman');