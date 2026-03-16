fn = '5m_range';
load(strcat([fn,'.mat']));
i = 1;
[~,max_position]=find( (data_s(i,:)>1e4)>0);
max_position=max_position(1);
% figure;plot(data_s(i,:));
sig=data_s(i:i+300,max_position:end);

index=1306;
r2=sig(:,index); 
smoothedData = movmean(r2, 60);
r2=double(r2)- smoothedData;
% r2=filter(RR_BPF,r2);
f1=sqrt(abs(MYIAA(single(r2) ) ));
f_real=(1:length(f1)/2)-1;f_real=f_real*10/length(f1);
hb=f_real*60;
index_hb=intersect(find(hb>=10),find(hb<=50));
index1=size(hb,2)-size(hb(hb>6),2);index2=size(hb(hb<40),2);
f_real=(1:length(f1)/2)-1;f_real=f_real*10/length(f1);
f1=(f1(1:length(f1)/2));f1=f1/max(f1(index1:index2));
    index_hb=intersect(find(hb>=6),find(hb<=40));
    [peak,locs]=findpeaks(f1(index_hb));
    [maxpeak,loc_peak]=max(peak);
    ind1=locs(loc_peak-1)+index_hb(1)-1;
    ind2=locs(loc_peak+1)+index_hb(1)-1;
    side=[f1(1:ind1) ;f1(ind2:end)];
    side_level=sum(side.^2)/size(side,1);
    side_level=sum(f1.^2)/size(f1,1);
    PIL5m=maxpeak^2/side_level;

load('nloshuxi.mat');
x=data_save(107,100:400);
phase_breath=angle(x);
phase_temp = [diff(unwrap( phase_breath)),0];
phase_unwrap = filter(RR_BPF,phase_temp);
f2=MYIAA(phase_unwrap.');f2=(f2(1:length(f2)/2));
light5m=r2/max(r2);
mmwave5m=phase_unwrap/pi;
figure;plot([1:301]/10,light5m,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 32]);
hold on;
plot([1:301]/10,mmwave5m,'-.','LineWidth',2,'Color',[.6, 0.5, .0]);
set(gcf, 'Position', [100, 100, 600, 250]);set(gca, 'FontSize', 12);
xlabel('Slow Time (s)');ylabel('Normalized waveform'); 
legend('Photonic','mmWave');
ylim([-1.3 1.3]);
% set(gcf, 'Position', [100, 100, 800, 500]);
f5m=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,f5m,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
hold on;
% plot(hb,20*log10(f2)-max(20*log10(f2)),'LineWidth',2);
set(gcf, 'Position', [100, 100, 600, 250]);set(gca, 'FontSize', 12);
xlim([6 60]);
ylim([-42 2])
hold on
real_rate_5m=17.25;
scatter(real_rate_5m,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
legend('Photonic','Real breath rate'); 
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');

fn = '10m_range';
load(strcat([fn,'.mat']));
i = 1;
[~,max_position]=find( (data_s(i,:)>1e4)>0);
max_position=max_position(1);
% figure;plot(data_s(i,:));
sig=data_s(i:i+300,max_position:end);

index=2650;
r2=sig(:,index); 
smoothedData = movmean(r2, 60);
r2=double(r2)- smoothedData;
f1=sqrt(abs(MYIAA(single(r2) ) ));
f_real=(1:length(f1)/2)-1;f_real=f_real*10/length(f1);

    index_hb=intersect(find(hb>=6),find(hb<=40));
    [peak,locs]=findpeaks(f1(index_hb));
    [maxpeak,loc_peak]=max(peak);
    ind1=locs(loc_peak-1)+index_hb(1)-1;
    ind2=locs(loc_peak+1)+index_hb(1)-1;
    side=[f1(1:ind1) ;f1(ind2:end)];
    side_level=sum(side.^2)/size(side,1);
    side_level=sum(f1.^2)/size(f1,1);
    PIL10m=maxpeak^2/side_level;
hb=f_real*60;
index_hb=intersect(find(hb>=10),find(hb<=50));
index1=size(hb,2)-size(hb(hb>6),2);index2=size(hb(hb<40),2);
f1=(f1(1:length(f1)/2));f1=f1/max(f1(index1:index2));
light_10m=r2/max(r2);
figure;plot([1:301]/10,light_10m,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 32]);ylim([-1.3 1.3]);
set(gcf, 'Position', [100, 100, 600, 250]);set(gca, 'FontSize', 12);
xlabel('Slow Time (s)');ylabel('Normalized waveform'); 
legend('Photonic');
% set(gcf, 'Position', [100, 100, 800, 500]);
f_10m=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,f_10m,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 60]);
ylim([-40 2])
hold on
real_rate_10=12.05;
scatter(real_rate_10,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
% set(gcf, 'Position', [100, 100, 800, 500]);
legend('Photonic', 'Real breath rate'); 
set(gcf, 'Position', [100, 100, 600, 250]);set(gca, 'FontSize', 12);


fn = 'handsmove2';
load(strcat([fn,'.mat']));
sig=data_p643;
index=643;
r2=sig(:,index); 
smoothedData = movmean(r2, 60);
r2=double(r2)- smoothedData;
f1=sqrt(abs(MYIAA((r2) ) ));
f_real=(1:length(f1)/2)-1;f_real=f_real*10/length(f1);
hb=f_real*60;
index_hb=intersect(find(hb>=10),find(hb<=50));
index1=size(hb,2)-size(hb(hb>6),2);index2=size(hb(hb<40),2);
f1=(f1(1:length(f1)/2));f1=f1/max(f1(index1:index2));
    index_hb=intersect(find(hb>=6),find(hb<=40));
    [peak,locs]=findpeaks(f1(index_hb));
    [maxpeak,loc_peak]=max(peak);
    ind1=locs(loc_peak-1)+index_hb(1)-1;
    ind2=locs(loc_peak+1)+index_hb(1)-1;
    side=[f1(1:ind1) ;f1(ind2:end)];
    side_level=sum(side.^2)/size(side,1);
    side_level=sum(f1.^2)/size(f1,1);
    PIL2m=maxpeak^2/side_level;
hb=f_real*60;
index_hb=intersect(find(hb>=10),find(hb<=50));
index1=size(hb,2)-size(hb(hb>6),2);index2=size(hb(hb<40),2);

load('radar_movehands.mat');
light_move=r2/max(r2);
figure;plot([1:301]/10,light_move,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 32]);ylim([-1.5 1.5]);
hold on;
a=1;
b=phase_unwrap(a:a+300);b=b/max(b);
mmwave_move=b;
plot([1:301]/10,mmwave_move,'-.','LineWidth',2,'Color', [.6, 0.5, .0])
legend('Photonic','mmWave');
xlabel('Slow Time (s)');ylabel('Normalized waveform');
set(gcf, 'Position', [100, 100, 600, 250]);set(gca, 'FontSize', 12);

f2=sqrt(abs(MYIAA(single(b') ) ));
f2=f2(1:length(f2)/2);f2=f2/max(f2);

f_move=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,f_move,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 60]);
ylim([-70 4]);
set(gca, 'FontSize', 12);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
hold on;
f_mmwave_move=20*log10(f2);
plot(hb,f_mmwave_move,'LineWidth',2.5);

xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
hold on
real_rate_move=22.2;
scatter(real_rate_move,0,40,'r^', 'LineWidth',2 ,'MarkerFaceColor', 'none');
legend('Photonic','mmWave', 'Real breath rate'); 
set(gcf, 'Position', [100, 100, 600, 250]);

%%

close all; 
set(0, 'DefaultAxesFontName', 'Arial'); % 
set(0, 'DefaultTextFontName', 'Arial');

% 
color_light = [0.12, 0.47, 0.71];   %  (Photonic)
color_mmwave = [0.84, 0.15, 0.16];  %  (mmWave)
color_real = [0.20, 0.63, 0.17];    % 
line_width = 1.2;
font_size_label = 10;
font_size_axis = 9;


set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultLineLineWidth', 1.2);

% 
color_photonic = [0, 0.4470, 0.7410]; % (Photonic)
color_gt = [0.8500, 0.3250, 0.0980];   %  (mmWave/GT)
color_spectrum = [0.4940, 0.1840, 0.5560]; %  (频谱线)
color_grid = [0.94, 0.94, 0.94];       % 

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 35, 15]);
tlo = tiledlayout(2, 3, 'TileSpacing', 'Loose', 'Padding', 'Compact');

% ---------------------------------------------------------

apply_style = @(ax, title_str) ...
    set(ax, 'Box', 'off', 'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', ...
    'GridColor', color_grid, 'GridAlpha', 1, 'FontSize',100, ...
    'LineWidth', 0.8, 'TitleFontSizeMultiplier', 1.1);

% ---------------------------------------------------------

ax1 = nexttile;
p1 = plot([1:301]/10, light5m, 'Color', color_photonic,'LineWidth',1.5); hold on;
p2 = plot([1:301]/10, mmwave5m, '--', 'Color', color_gt,'LineWidth',1.5);
ylabel('Norm. waveform'); xlabel('Time (s)');
% title('\textbf{a}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
legend([p1, p2], {'Photonic', 'mmWave'}, 'Box', 'off', 'Location', 'northeast');
apply_style(ax1); ylim([-1.3 1.3]); xlim([0 30]);


ax4 = nexttile(4);
plot(hb, f5m, 'Color', color_spectrum,'LineWidth',1.5); hold on;
scatter(real_rate_5m, 0, 45, color_gt, '^', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Norm. Amp. (dB)'); xlabel('Breaths per minute (bpm)');
% title('\textbf{d}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
legend('Photonic', 'Real rate', 'Box', 'off');
apply_style(ax4); xlim([6 60]); ylim([-42 5]);

ax2 = nexttile(2);
p1=plot([1:301]/10, light_10m, 'Color', color_photonic,'LineWidth',1.5);
ylabel('Norm. waveform'); xlabel('Time (s)');
legend([p1], {'Photonic'}, 'Box', 'off', 'Location', 'northeast');
% title('\textbf{b}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
apply_style(ax2); ylim([-1.3 1.3]); xlim([0 30]);


ax5 = nexttile(5);
plot(hb, f_10m, 'Color', color_spectrum,'LineWidth',1.5); hold on;
scatter(real_rate_10, 0, 45, color_gt, '^', 'Filled', 'MarkerEdgeColor', 'k');
legend('Photonic', 'Real rate', 'Box', 'off');
ylabel('Norm. Amp. (dB)'); xlabel('Breaths per minute (bpm)');
% title('\textbf{e}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
apply_style(ax5); xlim([6 60]); ylim([-40 5]);


ax3 = nexttile(3);
p1=plot([1:301]/10, light_move, 'Color', color_photonic,'LineWidth',1.5); hold on;
p2=plot([1:301]/10, mmwave_move, '--', 'Color', color_gt,'LineWidth',1.5);
ylabel('Norm. waveform'); xlabel('Time (s)');
% title('\textbf{c}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
legend([p1, p2], {'Photonic', 'mmWave'}, 'Box', 'off', 'Location', 'northeast');
apply_style(ax3); ylim([-1.5 1.5]); xlim([0 30]);


ax6 = nexttile(6);
plot(hb, f_move, 'Color', color_spectrum,'LineWidth',1.5); hold on;
plot(hb, f_mmwave_move, 'Color', color_gt, 'LineWidth', 1,'LineWidth',1.5);
scatter(real_rate_move, 0, 45, color_gt, '^', 'Filled', 'MarkerEdgeColor', 'k');
ylabel('Norm. Amp. (dB)'); xlabel('Breaths per minute (bpm)');
% title('\textbf{f}', 'Interpreter', 'latex', 'Units', 'normalized', 'Position', [-0.15, 1.05]);
legend('Photonic', 'mmWave', 'Real rate', 'Box', 'off', 'Location', 'northeast');
apply_style(ax6); xlim([6 60]); ylim([-75 10]);

set(findall(gcf,'-property','FontSize'),'FontSize',12);