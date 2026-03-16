data1=readtable('150cm_12_30bpm_f.csv');
v = table2array(data1(2:end,:));
v=v-mean(v);
v=v/max(v);
v=v(540:2600);
t=1/50:1/50:length(v)*1/50;
figure;
plot(t,v,'LineWidth',2.5);
load('150cm1.mat');
a1=a;


% --- 1. 
nodes = [
    0.02, 0.19, 0.38;  % 
    0.13, 0.40, 0.67;  % 
    0.40, 0.66, 0.81;  % 
    0.97, 0.97, 0.97;  % 
    0.94, 0.54, 0.38;  % 
    0.70, 0.09, 0.17;  % 
    0.40, 0.00, 0.12   % 
];
% 
custom_map = interp1(linspace(0, 1, size(nodes, 1)), nodes, linspace(0, 1, 256));

% --- 2.
figure('Color', 'w');
h = mesh(a); 

% 
shading interp; 

% --- 3. 
colormap(custom_map);
cb = colorbar;
ylabel(cb, 'Photon number', 'FontSize', 10, 'FontName', 'Arial');

% 
xl=xlabel('Range'); yl=ylabel('Slow time'); zlabel('Photon number');
set(gca, 'FontSize', 10, 'FontName', 'Arial', 'LineWidth', 1.1);
grid on; set(gca, 'GridAlpha', 0.1);
view(33, 30);
axis tight; % 
% 
set(gca, 'xtick', []); 
% set(xl, 'Rotation', -20);     
set(yl, 'Rotation', 38);      
% 
set(gca, 'ytick', []);
set(xl, 'VerticalAlignment', 'middle');
set(yl, 'VerticalAlignment', 'middle');

a=double(a( :,652));
smoothedData = movmean(a, 60);a=a-smoothedData ;
% a=a-mean(a);
a=a(110:end);
f1=sqrt(abs(MYIAA(single(a(100:end)) ) ));%
% f1=abs(fft(a(100:end),3070));
a=a/max(a);
t1=0.1:0.1:length(a)*0.1;
plot(t1,a*1.2,'LineWidth',2.5,'Color',[0, 0.4470, 0.7410]);hold on;
plot(t,v,'-.','LineWidth',2.5,'Color',[0.8500, 0.3250, 0.0980]);
set(gca, 'FontSize', 15);
xlim([0 45]);ylim([-2 1.8]);
xlabel('Slow Time (s)');ylabel('Normalized waveform');
legend('Light','GT');
set(gcf, 'Position', [100, 100, 1100, 300]);

for i=1:100 
a_resampled = resample(a, 5, 1);
v1=v(1:2030); 
R = corrcoef(v1,a_resampled );
RR(i)=R(1,2);
end


max_f=1/(100/1000)/2;%
f_real=(1:length(f1)/2)-1;
f_real=f_real*max_f/length(f_real); 
hb=f_real*60;
index_hb=intersect(find(hb>=6),find(hb<=40));
f1=(f1(1:length(f1)/2));
figure;plot(hb,20*log10(f1)-max(20*log10(f1)), 'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2.5);
xlim([6 60]);
ylim([-50 2]);
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
set(gcf, 'Position', [100, 100, 600, 350]);
hold on
real_rate=12.30;
p2=scatter(real_rate,0,40,'^', 'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2 ,'MarkerFaceColor', 'none');
legend(p2, 'Real breath rate'); 

%%


true_vals=[9.047 9.6463 11.4431 13.0760 13.7435 14.0093 17.2304 17.7515 19.9396 26.7346  27.3437 35.8040 ];
meas_vals=[8.9701 9.56811 11.7608 13.1561 14.3522 13.7542 17.3422 17.3422 19.5349 27.5083   27.907 36.0797 ];
n = numel(true_vals);
err = meas_vals - true_vals;
MAE = mean(abs(err));
RMSE = sqrt(mean(err.^2));
bias = mean(err);
[rho, pval] = corrcoef( meas_vals(:),true_vals(:)); 
if numel(rho)>1, pearson_r = rho(1,2); else pearson_r = NaN; end

% 
[~, idx] = sort(abs(meas_vals), 'ascend'); % 
t_sorted = true_vals(idx);
m_sorted = meas_vals(idx);

%
figure('Color','w','Position',[100 100 1200 300]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% 1) 
% nexttile([1 2]);
X = [t_sorted(:) m_sorted(:)];
hb2 = bar(X,'grouped','BarWidth',0.75);
hold on;
% 
ngroups = size(X,1); nbars = size(X,2);
for i = 1:nbars
    % XEndPoints 
    try
        xpos(i,:) = hb2(i).XEndPoints;
    catch
        % 
        x = (1:ngroups) - 0.15 + (2*(i-1))*(0.3/nbars); 
        xpos(i,:) = x;
    end
end
% 
for k = 1:ngroups
    plot(xpos(:,k), X(k,:), '-k', 'LineWidth', 1);
end
% 
legend({'Estimated','GT'}, 'Location','best');
xlabel('Trial (sorted by breath rate)');
ylabel('Breath per minute (bpm)');
% title(sprintf('Paired comparison (MAE=%.3g, RMSE=%.3g, bias=%.3g)', MAE, RMSE, bias));
set(gca,'FontSize',11,'Box','on','XTick',1:ngroups);

% 
for k = 1:ngroups
    xt = xpos(2,k); yt = X(k,2);
    text(xt, yt + 0.02*range([X(:)]), sprintf('%.2f', err(idx(k))), ...
        'HorizontalAlignment','center','FontSize',9);
end
hold off;
 

%%
 %%
close all;
% 
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 14]);

%
C_GT = [0.8500, 0.3250, 0.0980]; % GT 
C_Light = [0, 0.4470, 0.7410];  % Light 
C_IAA = [0.4940, 0.1840, 0.5560]; % IAA 

% 
left_col = 0.08; 
right_col = 0.65;
width_large = 0.88; % 
width_half = 0.5;  % 
%% --- Panel (c):
% 
ax_c = axes('Position', [left_col, 0.65,0.5, 0.28]);
t_v1 = (1:length(v))/50;
t_a = (1:length(a))*0.1;
p1 = plot(t_a, a*1.2, 'Color', C_Light, 'LineWidth', 1.8); hold on;
p2 = plot(t_v1, v, '-.', 'Color', C_GT, 'LineWidth', 1.5);
ylabel({'Norm.','Waveform'}); xlabel('Time (s)');
title('c', 'Units', 'normalized', 'Position', [-0.1, 1.05], 'FontSize', 14, 'FontWeight', 'bold');
legend([p1, p2], {'Photonic NLOS', 'GT'}, 'Box', 'off', 'FontSize', 9);
grid on; xlim([0 45]); ylim([-2 1.8]);
set(ax_c, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);

%% --- Panel (d): 
ax_d = axes('Position', [right_col, 0.65, 0.3, 0.28]);
plot(hb, 20*log10(f1)-max(20*log10(f1)), 'Color', C_IAA, 'LineWidth', 1.8); hold on;
color_gt = [0.8500, 0.3250, 0.0980]; 
p_real = scatter(real_rate, 0, 45,color_gt, '^', 'Filled','MarkerEdgeColor',  'k');
ylabel({'Norm.','Amp. (dB)'}); xlabel('Breaths per minute (bpm)');
title('d', 'Units', 'normalized', 'Position', [-0.15, 1.05], 'FontSize', 14, 'FontWeight', 'bold');
legend(p_real, 'Real breath rate', 'Box', 'off', 'FontSize', 9);
grid on; xlim([6 60]); ylim([-50 5]);
set(ax_d, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10);

%% --- Panel (e): 
ax_e = axes('Position', [left_col, 0.15, width_large, 0.35]);
X_data = [t_sorted(:) m_sorted(:)];
hb_bar = bar(X_data, 'grouped', 'BarWidth', 0.8, 'EdgeColor', 'none');
hb_bar(1).FaceColor = C_Light; % Estimated
hb_bar(2).FaceColor = C_GT;    % GT
hold on;

% 
for i = 1:numel(hb_bar)
    xpos_data(i,:) = hb_bar(i).XEndPoints;
end
for k = 1:size(X_data,1)
    plot(xpos_data(:,k), X_data(k,:), '-k', 'LineWidth', 0.5); % 
end

% 
for k = 1:size(X_data,1)
    text(xpos_data(2,k), X_data(k,2) + 1.5, sprintf('%.2f', err(idx(k))), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);
end

ylabel('Breaths per minute (bpm)'); xlabel('Trial (sorted by breath rate)');
title('e', 'Units', 'normalized', 'Position', [-0.06, 1.05], 'FontSize', 14, 'FontWeight', 'bold');
legend({'Estimated', 'GT'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 10);
set(ax_e, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 10, 'XTick', 1:size(X_data,1));
grid on; ylim([0 45]);

close all;
%% ---
close all;
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 28, 15]);

%
C_GT = [0.8500, 0.3250, 0.0980]; 
C_Light = [0, 0.4470, 0.7410];  
C_IAA = [0.4940, 0.1840, 0.5560]; 
top_row_y = 0.60;    % 
top_row_h = 0.30;    % 

%% --- Panel (a): 3D Mesh 
% Position: [Left, Bottom, Width, Height]
%% --- Panel (a) 
ax_a = axes('Position', [0.05, top_row_y, 0.252, top_row_h]);

% 
h3d = mesh(a1 / 1000); 

shading interp; colormap(ax_a, custom_map);
cb = colorbar;
cb.YAxisLocation = 'left'; %
% 
% ylabel(cb, 'Photon number (\times 10^3)', 'FontSize', 8, 'FontName', 'Arial');

view(28, 30); grid on; set(gca, 'GridAlpha', 0.1);
set(gca, 'xtick', [], 'ytick', [], 'FontSize', 9, 'FontName', 'Arial');

% 
xlabel('Range'); yl = ylabel('Slow time'); 
zl = zlabel('Photon number (\times 10^3)'); 

set(yl, 'Rotation', 42, 'VerticalAlignment', 'middle');
% title('d', 'Units', 'normalized', 'Position', [-0.12, 1.2], 'FontSize', 14, 'FontWeight', 'bold');
axis tight;

%% --- Panel (c): 
% 
ax_c = axes('Position', [0.332, top_row_y, 0.378, top_row_h]);
p1 = plot(t1, a*1.2, 'Color', C_Light, 'LineWidth', 1.5); hold on;
p2 = plot(t, v, '-.', 'Color', C_GT, 'LineWidth', 1.2);
ylabel({'Norm.','Waveform'}); xlabel('Time (s)');
% title('e', 'Units', 'normalized', 'Position', [-0.08, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
legend([p1, p2], {'Photonic NLOS', 'GT'}, 'Box', 'off', 'FontSize', 8, ...
    'Orientation', 'horizontal', 'Location', 'northoutside');
grid on; xlim([0 45]); ylim([-2 1.8]);
set(ax_c, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 9);

%% --- Panel (d): 
% 
ax_d = axes('Position', [0.74, top_row_y, 0.21, top_row_h]);
plot(hb, 20*log10(f1)-max(20*log10(f1)), 'Color', C_IAA, 'LineWidth', 1.5); hold on;
p_real = scatter(real_rate, 0, 35, C_GT, '^', 'Filled', 'LineWidth', 1.2,'MarkerEdgeColor',  'k');
ylabel({'Norm.','Amp. (dB)'}); xlabel('Breaths per minute (bpm)');
% title('f', 'Units', 'normalized', 'Position', [-0.18, 1.1], 'FontSize', 14, 'FontWeight', 'bold');
legend(p_real, 'Real rate', 'Box', 'off', 'FontSize', 8, 'Location', 'northeast');
grid on; xlim([6 60]); ylim([-50 5]);
set(ax_d, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 9);

%% --- Panel (e): 
% 
ax_e = axes('Position', [0.05, 0.12, 0.90, 0.35]);
hb_bar = bar([t_sorted(:) m_sorted(:)], 'grouped', 'BarWidth', 0.8, 'EdgeColor', 'none');
hb_bar(1).FaceColor = C_Light; hb_bar(2).FaceColor = C_GT;
hold on;
%
for i = 1:numel(hb_bar), xpos_data(i,:) = hb_bar(i).XEndPoints; end
for k = 1:size(t_sorted,2)
    plot(xpos_data(:,k), [t_sorted(k) m_sorted(k)], '-k', 'LineWidth', 0.5); 
    text(xpos_data(2,k), m_sorted(k) + 2, sprintf('%.2f', err(idx(k))), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontName', 'Arial');
end
ylabel('Breaths per minute (bpm)'); xlabel('Trial (sorted by breath rate)');
% title('g', 'Units', 'normalized', 'Position', [-0.03, 1.05], 'FontSize', 14, 'FontWeight', 'bold');
set(ax_e, 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 9, 'XTick', 1:numel(t_sorted));
grid on; ylim([0 48]);
legend({'Estimated', 'GT'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 10);