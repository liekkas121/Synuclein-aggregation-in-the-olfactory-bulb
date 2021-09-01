% Set default parameters
fs = 1000;
pre = 4;
post = 6;
stimu = 2;

% Select the folders to be analyzed
FolderList = uipickfiles;

%% Pre-analysis
for f = 1:size(FolderList,2)
    
    % Determine the current subfolder name
    selected_path = FolderList{f};
    cd(selected_path);
    work_folder = selected_path(find(selected_path == '\',1,'last')+1:end);
    disp([work_folder,' is being processed . . .']);
    
    % Determine mouse type
    unit_mice_type(f) = contains(work_folder,'Ctrl');
    
    % Load 'Odor¡¢Behavior' raw data from 'Behavior.txt'
    change_suffix('dat','txt');
    try
        odor = str2num(cell2mat(load_result('Behavior.txt',4)));
    catch
        odor = textread('Behavior.txt','%n');
    end
    
    % Find event
    event = find_event(2000);

    % Check number of Trials and Behaviors
    if length(event) ~= length(odor)
        error('Error>> The number of Trials is different from the number of Behaviors !')
    end
    
    % Delete event with noise
    if exist('Del.mat','file')~=0
        load('Del.mat');
        event(logical(LFP_del),:) = [];
        odor(logical(LFP_del),:) = [];
    end
    
    % Load lfp_raw
    lfp_index = 2;
    LFPnm = dir(['FP','*']);
    load(LFPnm(lfp_index).name);
    lfp_raw = eval(['FP',LFPnm(lfp_index).name(end-5:end-4)]);

    % Split 'lfp_raw' by 'event'
    for j = 1:length(event)
        lfp(j,:) = lfp_raw(event(j)-pre*fs+1:event(j)+post*fs);
    end
    
    % STFT
    bin = 0.1;
    window = hamming(1*fs);
    noverlap = (1-bin)*fs;
    nfft = max(256,2^nextpow2(length(window)));
    for k = 1:length(event)
        x = lfp(k,:);
        [~,freq,time,power(:,:,k)] = spectrogram(x,window,noverlap,nfft,fs);
    end
    
    % Normalization
    power = power./mean(mean(power(:,1:pre/bin,:),3),2);
    
    % Sort Odor
    % power_odor data structure : Mice*odor
    for odor_chan = min(odor):max(odor)
        power_odor{f,odor_chan+1} = power(freq < 100,:,odor == odor_chan);
    end

    clearvars -except fs pre post stimu FolderList f freq time unit_mice_type power_odor
end

%% Analysis
% Load mat file
load('E:\Documents\MATLAB\2020-07-02 OB-Syn-MT-Array Fengjiao Chen\02_Mat\LFP_Evoked_Synuclein.mat');

[theta,beta,loga,higa] = find_band(freq);

roi_mice = 1:18;

fig1 = figure('name','LFP_Evoked_Synuclein','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);

% Mean
odor_power_time(:,:,1) = cell2mat(cellfun(@(x) mean2(squeeze(mean(x(beta,time-pre>0 & time-pre <=2,:),1))),power_odor(roi_mice,:),'UniformOutput',false));
odor_power_time(:,:,2) = cell2mat(cellfun(@(x) mean2(squeeze(mean(x(higa,time-pre>0 & time-pre <=2,:),1))),power_odor(roi_mice,:),'UniformOutput',false));

y1(:,1) = reshape(odor_power_time( unit_mice_type(roi_mice),:,1),[],1);
y2(:,1) = reshape(odor_power_time(~unit_mice_type(roi_mice),:,1),[],1);
y1(:,2) = reshape(odor_power_time( unit_mice_type(roi_mice),:,2),[],1);
y2(:,2) = reshape(odor_power_time(~unit_mice_type(roi_mice),:,2),[],1);

% Bar_scatter of spotaneous_mfr
unpaired_stats1 = unpaired_test({y1(:,1),y2(:,1)});

subplot(3,3,1)
hold on
bar_scatter_alpha(1,y1(:,1),'r');
bar_scatter_alpha(2,y2(:,1),'b');
ylabel('Normalized Power','fontsize',10);
ylim([0 6]);
title({unpaired_stats1.method;[unpaired_stats1.stat,', ',unpaired_stats1.p]},'FontSize',6);

% High gamma
unpaired_stats2 = unpaired_test({y1(:,2),y2(:,2)});

subplot(3,3,2)
hold on
bar_scatter_alpha(1,y1(:,2),'r');
bar_scatter_alpha(2,y2(:,2),'b');
ylabel('Normalized Power','fontsize',10);
ylim([0 1.5]);
title({unpaired_stats2.method;[unpaired_stats2.stat,', ',unpaired_stats2.p]},'FontSize',6);

% Cumulative fraction units of spotaneous_mfr
[~,p,ks2stat] = kstest2(y1(:,1),y2(:,1));
subplot(3,3,4)
hold on
h1 = cdfplot(y1(:,1));
h2 = cdfplot(y2(:,1));
set( h1, 'LineStyle', '-', 'Color', 'r', 'linewidth', 0.75);
set( h2, 'LineStyle', '-', 'Color', 'b', 'linewidth', 0.75);
box off; title(''); grid off;
set(gca,'tickdir','out','ticklength',[0.025 0.025],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlim([0 6]);
xlabel('Normalized Power','fontsize',8);
ylabel('Cumulative fraction units','fontsize',8);
title({'Two-sample K-S test';[['Z', ' = ',num2str(ks2stat)],', ',['P = ',num2str(p)]]},'FontSize',6);

[~,p,ks2stat] = kstest2(y1(:,2),y2(:,2));
subplot(3,3,5)
hold on
h1 = cdfplot(y1(:,2));
h2 = cdfplot(y2(:,2));
set( h1, 'LineStyle', '-', 'Color', 'r', 'linewidth', 0.75);
set( h2, 'LineStyle', '-', 'Color', 'b', 'linewidth', 0.75);
box off; title(''); grid off;
set(gca,'tickdir','out','ticklength',[0.025 0.025],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlim([0 1.2]);
xlabel('Normalized Power','fontsize',8);
ylabel('Cumulative fraction units','fontsize',8);
title({'Two-sample K-S test';[['Z', ' = ',num2str(ks2stat)],', ',['P = ',num2str(p)]]},'FontSize',6);

%% Imagesc
odor_power_time_beta = cellfun(@(x) mean(mean(x(beta,:,:),1),3),power_odor(roi_mice,:),'UniformOutput',false);
odor_power_time_higa = cellfun(@(x) mean(mean(x(higa,:,:),1),3),power_odor(roi_mice,:),'UniformOutput',false);

y1_beta = cell2mat(reshape(odor_power_time_beta( unit_mice_type(roi_mice),:,1),[],1));
y2_beta = cell2mat(reshape(odor_power_time_beta(~unit_mice_type(roi_mice),:,1),[],1));
y1_higa = cell2mat(reshape(odor_power_time_higa( unit_mice_type(roi_mice),:,1),[],1));
y2_higa = cell2mat(reshape(odor_power_time_higa(~unit_mice_type(roi_mice),:,1),[],1));

subplot(341)
[~,y1_beta_I] = sort(y1(:,1));
imagesc(time-pre,1:size(y1_beta,1),y1_beta(y1_beta_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'fontsize',6);
set(gca,'ytick',[1,size(y1_beta,1)],'ytickLabel',{'1',size(y1_beta,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8); ylabel(('# of Animal-odor pairs'),'fontsize',8);
caxis([0,6]);

subplot(342)
[~,y2_beta_I] = sort(y2(:,1));
imagesc(time-pre,1:size(y2_beta,1),y2_beta(y2_beta_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'fontsize',6);
set(gca,'ytick',[1,size(y2_beta,1)],'ytickLabel',{'1',size(y2_beta,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8); ylabel(('# of Animal-odor pairs'),'fontsize',8);
caxis([0,6]);

subplot(343)
[~,y1_higa_I] = sort(y1(:,2));
imagesc(time-pre,1:size(y1_higa,1),y1_higa(y1_higa_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'fontsize',6);
set(gca,'ytick',[1,size(y1_higa,1)],'ytickLabel',{'1',size(y1_higa,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8); ylabel(('# of Animal-odor pairs'),'fontsize',8);
caxis([0,1.5]);

subplot(344)
[~,y2_higa_I] = sort(y2(:,2));
imagesc(time-pre,1:size(y2_higa,1),y2_higa(y2_higa_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'fontsize',6);
set(gca,'ytick',[1,size(y2_higa,1)],'ytickLabel',{'1',size(y2_higa,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8); ylabel(('# of Animal-odor pairs'),'fontsize',8);
caxis([0,1.5]);



