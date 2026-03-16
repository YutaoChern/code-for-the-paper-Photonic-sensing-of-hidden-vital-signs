clc;clear all;close all; %% 
T=33;T1=T;
tic
data_sum=zeros(16,16*T,4000);
for i=1:T
fn = 'data0926_dingwei\data';
load([fn,num2str(i),'.mat']);
data=zeros(16,16,4000);
for k2=1:16
        if mod(k2,2)~=0
            data(1:16,k2,:)=output([1:16]+16*(k2-1),:); %one point 4ms
        else
            data(16:-1:1,k2,:)=output([1:16]+16*(k2-1),:);
        end
end
data_sum(:,(i-1)*16+1:16*i,:)= data;
end
data_reshaped = reshape(data_sum, 16, 4,4*T, 4000);  % 
data = squeeze(sum(data_reshaped, 2));  % 

% plot(squeeze(data(1,1,:)));
sig=data;
blur = reshape([1:3,3:-1:1],[6,1,1]);
tempsig = convn(sig,blur,'same');
newsig=[];
for i = 1:16
            if i == 1
                    [peaks,loc1] =findpeaks(squeeze(data(1,1,1:1000)),'MinPeakHeight',100);
                    figure;plot(squeeze(data(1,1,:)));
                    pic1(i) =loc1(1);
            else
                [peaks,loc1] =findpeaks(squeeze(data(i,1,1:1000)),'MinPeakHeight',100);
                pic1(i) = loc1(1);
            end
            newsig(i,:,:)=sig(i,:,pic1(i):pic1(i)+2300);
end
newsig(newsig<0)=0;
% figure;plot(squeeze(newsig(10,10,:)));
g1 = 150;
newsig(:,:,1:g1) = 0;
figure;plot(squeeze(sum(sum(abs(newsig),1),2)));
data = zeros(size(data));
data(:,:,1:2301)=newsig;
data_target=data;

%%
T=20;
data_sum=zeros(16,16*T,4000);
for i=1:T
fn = 'data0926_kong\data';
load([fn,num2str(i),'.mat']);
data=zeros(16,16,4000);
for k2=1:16
    for k1=1:16
        if mod(k2,2)~=0
            data(k1,k2,:)=output(k1+16*(k2-1),:); %
        else
            data(k1,k2,:)=output(16-k1+1+16*(k2-1),:);
        end
    end
end
data_sum(:,(i-1)*16+1:16*i,:)= data;
end
data_reshaped = reshape(data_sum, 16, 4, T*4, 4000);  % 
data = squeeze(sum(data_reshaped, 2));  % 

sig=data;
blur = reshape([1:3,3:-1:1],[6,1,1]);
tempsig = convn(sig,blur,'same');
newsig=[];
for i = 1:16
            if i == 1
                    [peaks,loc1] =findpeaks(squeeze(data(1,1,1:1000)),'MinPeakHeight',100);
                    pic1(i) =loc1(1);
            else
                %ind = ind+1;
                [peaks,loc1] =findpeaks(squeeze(data(i,1,1:1000)),'MinPeakHeight',100);
                pic1(i) = loc1(1);
            end
            newsig(i,:,:)=sig(i,:,pic1(i):pic1(i)+2300);
end
newsig(newsig<0)=0;
figure;plot(squeeze(newsig(10,10,:)));
g1 = 150;
newsig(:,:,1:g1) = 0;
figure;plot(squeeze(sum(sum(abs(newsig),1),2)));
data = zeros(size(data));
data(:,:,1:2301)=newsig;
data=(sum(data,2)/80);
data_no=data;
%%
data1=data_target-data_no;
data1(data1<0)=0;
figure;plot(squeeze(sum(sum(abs(data1),1),2)));

% for i=1:T1*4
% data2=squeeze(data1(:,i,:));
% % figure;plot(squeeze((sum(abs(data2),1))));
% % figure;plot(movmean(squeeze((sum(abs(data2),1))),50));
% % figure;findpeaks(movmean(squeeze((sum(abs(data2),1))),50),'MinPeakWidth',30,'MinPeakHeight',50) ;
% [peks,locs]=findpeaks(movmean(squeeze((sum(abs(data2),1))),50),'MinPeakWidth',30,'MinPeakHeight',50) ;
% locs_T(i,:)=locs;
% end


binwidth=20*3/2*1e-4;
[x0] = linspace(-1,1,16);
y0=0;
nx=256;nz=1000;
rho_T=zeros(4*T1,nx,nz);
for i=1:4*T1
    rho = 1;
    data2=data1(:,i,:);
for j = 1:16
        data_temp = squeeze(data2(j,1,:));
        s = scale_my(size(data_temp,1),1000,1000);
        data_temp = data_temp.*s; % range 
        % plot(data_temp);pause(0.2) % 
        xstart=-3;xend=3;zstart=0;zend=6;
        [rho_temp,X,Z] = rho_bp_2d_optimized(data_temp,binwidth,x0(j),y0,0,nx,1,nz,xstart,0,zstart,xend,0,zend);
        rho = rho+rho_temp;
end
rho_T(i,:,:)=rho;
end
%%
nx=256;nz=1000;xstart=-3;xend=3;zstart=0;zend=6;
figure;mesh(squeeze(rho_T(1,:,:)));
peakx=[95 197 200 139];
peakz=[224 350 491 620 ];
figure;
set(gcf, 'Color', 'k', 'InvertHardcopy', 'off'); % 
ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w'); % 
hold(ax, 'on');
scatter(peakx/nx*(xend-xstart)+xstart,peakz/nz*zend,200, 'rx', 'LineWidth', 2); %
xlim([-2.5 (2.5)]);
ylim([zstart (zend)-1]);
set(ax, 'GridColor', [0.5 0.5 0.5], 'FontSize', 12); 
grid on;
box on;
xlabel('X Position (m)', 'Color', 'w', 'FontSize', 14);
ylabel('Z Position (m)', 'Color', 'w', 'FontSize', 14);
title('Localization', 'Color', 'w', 'FontSize', 16);


for t=2:4*T1
    rho_temp=squeeze(rho_T(t,:,:));
    xfanwei=7;zfanwei=10;
for i=1:4
    mask=zeros(size(rho_temp));
    mask(peakx(t-1,i)-xfanwei:peakx(t-1,i)+xfanwei,peakz(t-1,i)-zfanwei:peakz(t-1,i)+zfanwei)=1;
    mask_data=mask.*rho_temp;
    max_value=max(mask_data(:));
    [peakx(t,i),peakz(t,i)]=find(mask_data==max_value,1, 'first');
  % peakvalue=smoothed_data(peakx(i),peakz(i));
end
end
toc
%%

peakx_T_smoothed = movmean(peakx,10, 1)/nx*6-3;
peakz_T_smoothed = movmean(peakz,10, 1)/nz*6;
dis_r=sqrt((peakx_T_smoothed(:,1)-peakx_T_smoothed(:,2)).^2+(peakz_T_smoothed(:,1)-peakz_T_smoothed(:,2)).^2);
figure('Color', 'k');
ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w');
hold(ax, 'on');

traj1 = plot(ax, NaN, NaN, 'r-', 'LineWidth', 1.5, 'Color', [1, 0.5, 0.5]); % 1
traj2 = plot(ax, NaN, NaN, 'b-', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 1]); % 2
traj3 = plot(ax, NaN, NaN, 'g-', 'LineWidth', 1.5, 'Color', [0.5, 1, 0.5]); % 3
traj4 = plot(ax, NaN, NaN, 'm-', 'LineWidth', 1.5, 'Color', [1, 0.5, 1]);   % 4

h1 = plot(ax, NaN, NaN, 'ro', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerFaceColor', 'r');
h2 = plot(ax, NaN, NaN, 'bo', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
h3 = plot(ax, NaN, NaN, 'go', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerFaceColor', 'g');
h4 = plot(ax, NaN, NaN, 'mo', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerFaceColor', 'm');

xlabel('x (m)'); ylabel('z (m)'); title('Real-time Trajectory', 'Color', 'w');
legend({'TrajA', 'TrajB', 'TrajC', 'TrajD', 'A', 'B', 'C', 'D'}, ...
       'TextColor', 'w', 'Color', 'k', 'Location', 'best');
grid on;
xlim([-2.5, 2.5]); ylim([0,5]);

traj1_x = []; traj1_y = [];
traj2_x = []; traj2_y = [];
traj3_x = []; traj3_y = [];
traj4_x = []; traj4_y = [];

gif_filename = 'trajectory_animation.gif';
frame_rate = 100;
delay_time = 1/frame_rate; 

for k = 1:4*T1
    set(h1, 'XData', peakx_T_smoothed(k,1), 'YData', peakz_T_smoothed(k,1));
    set(h2, 'XData', peakx_T_smoothed(k,2), 'YData', peakz_T_smoothed(k,2));
    set(h3, 'XData', peakx_T_smoothed(k,3), 'YData', peakz_T_smoothed(k,3));
    set(h4, 'XData', peakx_T_smoothed(k,4), 'YData', peakz_T_smoothed(k,4));
    
    traj1_x = [traj1_x; peakx_T_smoothed(k,1)];
    traj1_y = [traj1_y; peakz_T_smoothed(k,1)];
    set(traj1, 'XData', traj1_x, 'YData', traj1_y);
    
    traj2_x = [traj2_x; peakx_T_smoothed(k,2)];
    traj2_y = [traj2_y; peakz_T_smoothed(k,2)];
    set(traj2, 'XData', traj2_x, 'YData', traj2_y);
    
    traj3_x = [traj3_x; peakx_T_smoothed(k,3)];
    traj3_y = [traj3_y; peakz_T_smoothed(k,3)];
    set(traj3, 'XData', traj3_x, 'YData', traj3_y);
    
    traj4_x = [traj4_x; peakx_T_smoothed(k,4)];
    traj4_y = [traj4_y; peakz_T_smoothed(k,4)];
    set(traj4, 'XData', traj4_x, 'YData', traj4_y);
    
    title(sprintf('Time = %.1f s', k*16 * 16/1000), 'Color', 'w');
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    
    if k == 1
        imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
    
    drawnow;
end
title('Localization & Tracking', 'Color', 'w', 'FontSize', 16);
set(gca, 'Color', 'k'); 
set(gcf, 'InvertHardcopy', 'off');  % 
saveas(gcf, 'black_plot.png');
%%
% 
%  target1=[-0.18 3.7-> 1.06 3.67];
%  target2=[0.947 3.012 -> -0.933 2.68 ->  -1.066 1.87 ]
% target3=[ 1.439 2.244 -> -0.33 1.947 -> -0.35 1.1478 ]
% target4=[ -0.989 1.6236 -> 1.868 1.617];
figure;
set(gcf, 'Color', 'k', 'InvertHardcopy', 'off'); % 
ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w'); % 
hold(ax, 'on');
% 绘制红色叉号（线宽加粗）
scatter(-0.18,3.7,200, 'rx', 'LineWidth', 2); % 
scatter(0.947,3.012,200, 'rx', 'LineWidth', 2); %
scatter(1.439,2.244,200, 'rx', 'LineWidth', 2); % 
scatter(-0.989,1.6236 ,200, 'rx', 'LineWidth', 2); %
h=quiver(-0.18, 3.7, 1.06+0.18, 3.67-3.7, ...
       0, ...                  
       'LineStyle', '--', ... 
       'MaxHeadSize', 0.5, ... 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h.Head.LineStyle = 'solid';  % 
h.Head.LineWidth = 1.5;      % 
h=quiver(0.947, 3.012,  -0.933-0.947,2.88-3.012, ...
       0, ...                  % 
       'LineStyle', '--', ...  % 
       'MaxHeadSize', 0.0, ... % 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h=quiver(-0.933, 2.88,  -1.066+0.933,1.87-2.88, ...
       0, ...                  % 
       'LineStyle', '--', ...  % 
       'MaxHeadSize', 0.5, ... % 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h.Head.LineStyle = 'solid';  % 
h.Head.LineWidth = 1.5;      % 
h=quiver(1.439, 2.244,   -0.33-1.439, 1.997-2.244, ...
       0, ...                  % 
       'LineStyle', '--', ...  % 
       'MaxHeadSize', 0.0, ... % 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h=quiver(-0.33, 1.947,  -0.35+0.33, 1.1478-1.997 , ...
       0, ...                  % 
       'LineStyle', '--', ...  % 
       'MaxHeadSize', 0.5, ... % 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h.Head.LineStyle = 'solid';  % 
h.Head.LineWidth = 1.5;      % 
h=quiver(-0.989, 1.6236, 1.868+0.989, 1.617-1.6236, ...
       0, ...                  % 
       'LineStyle', '--', ...  % 
       'MaxHeadSize', 0.2, ... % 
       'Color', 'w', ...       % 
       'LineWidth', 1.5);      % 
h.Head.LineStyle = 'solid';  % 
h.Head.LineWidth = 1.5;      % 
xlim([-2.5 (2.5)]);
ylim([zstart (zend)-1]);
set(ax, 'GridColor', [0.5 0.5 0.5], 'FontSize', 12); % 
grid on;
box on;
xlabel('X Position (m)', 'Color', 'w', 'FontSize', 14);
ylabel('Z Position (m)', 'Color', 'w', 'FontSize', 14);
title('Localization', 'Color', 'w', 'FontSize', 16);


data=flip(permute(rho_T,[1 3 2]),2);
figure;
colormap(jet); % 
load("color.mat");
colormap(CustomColormap);
h = imagesc(squeeze(data(1,:,:)));
axis tight manual;
title('t: 1');

filename = 'animation.gif';
delayTime = 0.05; 

for t = 1:size(data,1)
    set(h, 'CData', squeeze(data(t,:,:)));
    title(['t: ' num2str(t)]);
    drawnow;
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    
    if t == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end
%%

colors = [
    0.00, 0.45, 0.74;  % A: 
    0.85, 0.33, 0.10;  % B: 
    0.47, 0.67, 0.19;  % C: 
    0.49, 0.18, 0.56   % D: 
];


fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 20, 10]);
ax = axes('Box', 'off', 'TickDir', 'out', 'LineWidth', 1.0, 'FontName', 'Arial');
hold on; grid on;
set(ax, 'GridAlpha', 0.1, 'GridLineStyle', '--');

h_traj = gobjects(4, 1);   %
h_target = gobjects(4, 1); % 

for i = 1:4
    h_traj(i) = plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5);

    h_target(i) = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
                       'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
end


h_text = gobjects(4, 1);
tags = {'A', 'B', 'C', 'D'};
for i = 1:4
    h_text(i) = text(NaN, NaN, tags{i}, 'FontSize', 12 , ...
                     'Color', colors(i,:), 'FontName', 'Arial');
end

legend_handles = [h_traj; h_target];
legend_labels = {'Traj. A', 'Traj. B', 'Traj. C', 'Traj. D', ...
                 'Target A', 'Target B', 'Target C', 'Target D'};

lgd = legend(legend_handles, legend_labels, ...
             'NumColumns', 2, ...       % 
             'Box', 'off', ...          % 
             'Location', 'northeastoutside', ...
             'FontSize', 12, ...
             'FontName', 'Arial');

xlabel('X direction (m)' , 'FontSize', 14);
ylabel('Z direction (m)' , 'FontSize', 14);
xlim([-2.5, 2.5]); ylim([0, 5]);
daspect([1  1.7 1]); %

traj_x = cell(4,1); traj_z = cell(4,1);
for k = 1:4*T1
    for i = 1:4
        cx = peakx_T_smoothed(k, i);
        cz = peakz_T_smoothed(k, i);
        
        traj_x{i} = [traj_x{i}; cx];
        traj_z{i} = [traj_z{i}; cz];
        
        set(h_traj(i), 'XData', traj_x{i}, 'YData', traj_z{i});
        set(h_target(i), 'XData', cx, 'YData', cz);
        set(h_text(i), 'Position', [cx + 0.12, cz + 0.12]); %
    end
    
    title(sprintf('Localization & Tracking' ), 'FontSize', 14);
    drawnow;
end