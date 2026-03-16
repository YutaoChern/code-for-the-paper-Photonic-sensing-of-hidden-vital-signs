
%%
load('data0926_four.mat');
d1=max(data(150:450,935:944).');
d2=max(data(150:450,1205:1209).');
d3=max(data(200:500,1385:1395).');
d4=data(50:350,1542).';
subplot(4,1,1);plot(d1);
subplot(4,1,2);plot(d2);
subplot(4,1,3);plot(d3);
subplot(4,1,4);plot(d4);

smoothedData = movmean(d1, 60);
r2=d1-smoothedData;
ra=r2;
f1=sqrt(abs(MYIAA(single(r2.') ) ));%
max_f=1/(16*8/1000)/2;% smaple iterval 16*8 ms，corresponding max frequ=1/(delta_t)/2
f_real=(1:length(f1)/2)-1;f_real=f_real*max_f/length(f_real); %transfer to real fre
hb=f_real*60;
index_hb=intersect(find(hb>=6),find(hb<=40));
f1=(f1(1:length(f1)/2));
fa1=f1;
fa=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,20*log10(f1)-max(20*log10(f1)),'Color',  [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 120]);
ylim([-50 0]);
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
set(gcf, 'Position', [100, 100, 750, 350]);
% hold on
% real_rate=24.50;
% p2=scatter(real_rate,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
% legend(p2, 'Real respiratory rate'); 
figure;plot([1:301]*0.128,r2,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 40]);
set(gca, 'FontSize', 15);
xlabel('Time (s)');ylabel('Photon number');
set(gcf, 'Position', [100, 100, 750, 350]);

smoothedData = movmean(d2, 60);
r2=d2-smoothedData;
rb=r2;
f1=sqrt(abs(MYIAA(single(r2.') ) ));%
max_f=1/(16*8/1000)/2;% 
f_real=(1:length(f1)/2)-1;f_real=f_real*max_f/length(f_real); %
hb=f_real*60;
index_hb=intersect(find(hb>=6),find(hb<=40));
f1=(f1(1:length(f1)/2));
fb1=f1;
fb=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,20*log10(f1)-max(20*log10(f1)),'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 120]);
ylim([-50 0]);
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
set(gcf, 'Position', [100, 100, 750, 350]);
% hold on
% real_rate=31.40;
% p2=scatter(real_rate,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
% legend(p2, 'Real respiratory rate'); 
figure;plot([1:301]*0.128,r2,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 40]);
set(gca, 'FontSize', 15);
xlabel('Time (s)');ylabel('Photon number');
set(gcf, 'Position', [100, 100, 750, 350]);

smoothedData = movmean(d3, 60);
r2=d3-smoothedData;
rc=r2;
f1=sqrt(abs(MYIAA(single(r2.') ) ));%
max_f=1/(16*8/1000)/2;% 
f_real=(1:length(f1)/2)-1;f_real=f_real*max_f/length(f_real); %
hb=f_real*60;
index_hb=intersect(find(hb>=6),find(hb<=40));
f1=(f1(1:length(f1)/2));
fc1=f1;
fc=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,20*log10(f1)-max(20*log10(f1)),'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 120]);
ylim([-50 0]);
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
set(gcf, 'Position', [100, 100, 750, 350]);
% hold on
% real_rate=16.85;
% p2=scatter(real_rate,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
% legend(p2, 'Real respiratory rate'); 
figure;plot([1:301]*0.128,r2,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 40]);
set(gca, 'FontSize', 15);
xlabel('Time (s)');ylabel('Photon number');
set(gcf, 'Position', [100, 100, 750, 350]);

smoothedData = movmean(d4, 60);
r2=d4-smoothedData;
rd=r2;
f1=sqrt(abs(MYIAA(single(r2.') ) ));%
max_f=1/(16*8/1000)/2;% 
f_real=(1:length(f1)/2)-1;f_real=f_real*max_f/length(f_real); 
hb=f_real*60;
index_hb=intersect(find(hb>=6),find(hb<=40));
f1=(f1(1:length(f1)/2));
fd1=f1;
fd=20*log10(f1)-max(20*log10(f1));
figure;plot(hb,20*log10(f1)-max(20*log10(f1)),'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
xlim([6 120]);
ylim([-50 0]);
set(gca, 'FontSize', 15);
xlabel('Breaths per minute (bpm)');ylabel('Normalized amplitude (dB)');
set(gcf, 'Position', [100, 100, 750, 350]);
% hold on
% real_rate=26.78;
% p2=scatter(real_rate,0,40,'r','^', 'LineWidth', 1.2 ,'MarkerFaceColor', 'none');
% legend(p2, 'Real respiratory rate'); 
figure;plot([1:301]*0.128,r2,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);xlim([0 40]);
set(gca, 'FontSize', 15);
xlabel('Time (s)');ylabel('Photon number');
set(gcf, 'Position', [100, 100, 750, 350]);

%figure;plot(fa);hold on;plot(fb);hold on;plot(fc);hold on ;plot(fd);xlim([1 400]);
figure;
subplot(4,1,1);plot([1:301]*0.128,ra,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);ylim([-300 300]);
legend('Target A');
subplot(4,1,2);plot([1:301]*0.128,rb,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);ylim([-200 200]);
legend('Target B');
subplot(4,1,3);plot([1:301]*0.128,rc,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
legend('Target C');
subplot(4,1,4);plot([1:301]*0.128,rd,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);ylim([-100 100]);
legend('Target D');
han = axes('visible', 'off'); % 
han.XLabel.Visible = 'on';    %
han.YLabel.Visible = 'on';    % 
xlabel(han, 'Time (s)', 'FontSize', 15 );
ylabel(han, 'Photon number', 'FontSize', 15 );
set(gcf, 'Position', [00, 00,400, 700]);
%%

figure;
subplot(4,1,1);plot(hb,fa,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
[maxValue, maxIndex] = max(fa);
maxX = fa(maxIndex);
hold on; plot(hb(maxIndex), maxValue, 'o', 'MarkerSize', 5,'LineWidth',2);
text(hb(maxIndex) + 2.5,  maxValue-0.5, sprintf('%.2f', hb(maxIndex)), ...
    'VerticalAlignment', 'middle',  'FontWeight', 'bold');
xlim([6 120]);ylim([-40 0]);legend('Target A');
subplot(4,1,2);plot(hb,fb,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
[maxValue, maxIndex] = max(fb);
maxX = fa(maxIndex);
hold on; plot(hb(maxIndex), maxValue, 'o', 'MarkerSize', 5,'LineWidth',2);
text(hb(maxIndex) + 2.5,  maxValue-0.5, sprintf('%.2f', hb(maxIndex)), ...
    'VerticalAlignment', 'middle',  'FontWeight', 'bold');
xlim([6 120]);ylim([-40 0]);legend('Target B');
subplot(4,1,3);plot(hb,fc,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
[maxValue, maxIndex] = max(fc);
maxX = fa(maxIndex);
hold on; plot(hb(maxIndex), maxValue, 'o', 'MarkerSize', 5,'LineWidth',2);
text(hb(maxIndex) + 2.5,  maxValue-0.5, sprintf('%.2f', hb(maxIndex)), ...
    'VerticalAlignment', 'middle',  'FontWeight', 'bold');
xlim([6 120]);ylim([-40 0]);legend('Target C');
subplot(4,1,4);plot(hb,fd,'Color', [0.4940, 0.1840, 0.5560],'LineWidth',2);
[maxValue, maxIndex] = max(fd(1:500));
maxX = fa(maxIndex);
hold on; plot(hb(maxIndex), maxValue, 'o', 'MarkerSize', 5,'LineWidth',2);
text(hb(maxIndex) +2.5,  maxValue-0.5, sprintf('%.2f', hb(maxIndex)), ...
    'VerticalAlignment', 'middle',  'FontWeight', 'bold');
xlim([6 120]);ylim([-40 0]);legend('Target D');
han = axes('visible', 'off'); % 
han.XLabel.Visible = 'on';    % 
han.YLabel.Visible = 'on';    % 
xlabel(han, 'Breaths per minute (bpm)', 'FontSize', 15 );
ylabel(han, 'Normalized amplitude (dB)', 'FontSize', 15 );
set(gcf, 'Position', [00, 00,400, 700]);


hb_rand=intersect(find(6<hb),find(hb<40));
max(fd1(hb_rand))
A1= (max(fa1(hb_rand)).^2/mean(fa1.^2));
A2=max(fb1(hb_rand)).^2/mean(fb1.^2);
A3=max(fc1(hb_rand)).^2/mean(fc1.^2);
A4= (max(fd1(hb_rand)).^2/mean(fd1.^2));

%%
%% --- Figure 1: Waveforms (Time Domain) ---
%% --- Figure 1: Waveforms (Refined Y-Axis) ---
data_time = {ra, rb, rc, rd};
titles = {'Target A', 'Target B', 'Target C', 'Target D'};

yticks_cell = {[-300, 0, 300], [-200, 0, 200], [-100, 0, 100], [-100, 0, 100]};
ylims = {[-350 350], [-250 250], [-150 150], [-150 150]}; % 

fig1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 10, 18]);

tlo1 = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:4
    ax = nexttile;
    plot([1:301]*0.128, data_time{i}, 'LineWidth', 1.5, 'Color', [0.12, 0.47, 0.71]);
    hold on; grid on;
    set(ax, 'Box', 'off', 'TickDir', 'out', 'GridAlpha', 0.1, 'FontSize', 9);
    ylim(ylims{i});
    % 
    set(ax, 'YTick', yticks_cell{i}); 
    
    % 
    xlim([0, 301*0.128]);
    
    text(0.98, 0.82, titles{i}, 'Units', 'normalized', ...
        'FontSize', 12, 'HorizontalAlignment', 'right');
    
    if i < 4
        set(ax, 'XTickLabel', []);
    else
        xlabel('Time (s)', 'FontSize', 12 );
    end
end

ylabel(tlo1, 'Photon number', 'FontSize', 12, 'FontName', 'Arial' );

data_freq = {fa, fb, fc, fd};
fig2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [12, 2, 10, 18]);
tlo2 = tiledlayout(4, 1, 'TileSpacing', 'none', 'Padding', 'compact');

color_spec = [0.44, 0.19, 0.63]; % 
color_peak = [0.85, 0.33, 0.10]; % 
ylims = {[-45 5], [-45 5], [-45 5], [-35 15]}; % 
for i = 1:4
    ax = nexttile;
    curr_f = data_freq{i};
    
    % 
    plot(hb, curr_f, 'Color', color_spec, 'LineWidth', 1.2);
    hold on;
    
    % 
    range_idx = find(hb >= 6 & hb <= 110);
    [maxValue, localIdx] = max(curr_f(range_idx));
    peakX = hb(range_idx(localIdx));
    
    % 
    plot(peakX, maxValue, 'o', 'MarkerSize', 4, 'MarkerFaceColor', color_peak, 'MarkerEdgeColor', 'none');
    
    text(peakX + 3, maxValue + 1, sprintf('%.1f bpm', peakX), ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', color_peak);
    
    % 
    grid on;
    set(ax, 'Box', 'off', 'TickDir', 'out', 'GridAlpha', 0.1, 'FontSize', 9);
    xlim([6 120]);ylim(ylims{i});
    
    % 
    text(0.98, 0.82, titles{i}, 'Units', 'normalized', ...
         'FontSize', 12, 'HorizontalAlignment', 'right');
    
    if i < 4
        set(ax, 'XTickLabel', []);
    end
end

% 
xlabel(tlo2, 'Breaths per minute (bpm)', 'FontSize', 12 );
ylabel(tlo2, 'Norm. Amp. (dB)', 'FontSize', 12 );