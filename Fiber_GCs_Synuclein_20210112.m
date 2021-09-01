% Set default parameters
fs = 100;
pre = 5;
post = 10;
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
    
    % Mouse type
    con_mice(f) = contains(work_folder,'con');

    % Load 'Event、Lick、Fiber' raw data
    [ConvertedData,~,~,~,~]=convertTDMS(0,'Fiber.tdms');
    event_raw = ConvertedData.Data.MeasuredData(3).Data;
    lick_raw  = ConvertedData.Data.MeasuredData(5).Data;
    fiber_raw = ConvertedData.Data.MeasuredData(7).Data;
    
    % Load 'Odor、Behavior' raw data from 'Behavior.txt'
    odor = str2num(cell2mat(load_result('Behavior.dat',4)));
    
    % Find event
    event = find(diff(event_raw) > 0.9);
    
    % Debug
    if contains(work_folder,'syn13')
        event = event(2:end);
    end
        
    % Check number of Trials and Behaviors
    if length(event) ~= length(odor)
        error('Error>> The number of Trials is different from the number of Behaviors !')
    end
        
    % Split 'fiber_raw' by 'event'
    for j = 1:length(event)
        fiber(j,:) = fiber_raw(event(j)-pre*fs+1:event(j)+post*fs);
        lick(j,:) = lick_raw(event(j)-pre*fs+1:event(j)+post*fs);
    end
    
    % Normalization
    fiber = (fiber-mean(fiber(:,1:pre*fs),2))./mean(fiber(:,1:pre*fs),2)*100;
    
    % Sort Odor
    for odor_chan = min(odor):max(odor)
        animal_odor_fiber{f,odor_chan+1} = fiber(odor == odor_chan,:);
    end

        
    clearvars -except fs pre post bin stimu FolderList f con_mice animal_odor_fiber
end

%% Plot and imagesc single mouse         
time = -pre+1/fs : 1/fs : post;

fig1 = figure('name','Fiber_GCs_Synuclein_Imagesc_Plot_Shadow','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);

odor_fiber = animal_odor_fiber(1,:);

% Imagesc
for odor_chan = 1:size(odor_fiber,2)
    subplot(4,4,odor_chan);
    imagesc(time,1:size(odor_fiber{odor_chan},1),odor_fiber{odor_chan});
    colormap hot;
    
    set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
    xlabel('Time (s)','fontsize',8);ylabel(('Trial #'),'fontsize',6);
    set(gca,'xtick',[-4.99 0 2 5 10],'xtickLabel',{'-5','0','2','5','10'},'fontsize',6);
    set(gca,'ytick',[1,size(odor_fiber{odor_chan},1)],'ytickLabel',{'1',num2str(size(odor_fiber{odor_chan},1))},'fontsize',6);
    %caxis([-0.2,1]);
end

% Plot shadow
for odor_chan = 1:size(odor_fiber,2)
    subplot(4,4,odor_chan+8);
    plot_shadow(time,mean(odor_fiber{odor_chan},1),sem(odor_fiber{odor_chan},1),'k');
    set(gca,'tickdir','out','ticklength',[0.03 0.03],...
        'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
    xlabel('Time (s)','fontsize',10); ylabel(('ΔF/F (%)'),'fontsize',6);
    set(gca,'xtick',[-4.99 0 2 5 10],'xtickLabel',{'-5','0','2','5','10'},'fontsize',6);
    %ylim([-0.2,1]);
end

%% Statistic and plot all mice -- Mean
fig2 = figure('name','Fiber_GCs_Synuclein_Fiber_Statistic_Plot','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);

time = -pre+1/fs : 1/fs : post;

animal_odor_fiber_time = cellfun(@(x) mean(x,1),animal_odor_fiber,'UniformOutput',false);

y1_time = cell2mat(reshape(animal_odor_fiber_time( con_mice',:),[],1));
y2_time = cell2mat(reshape(animal_odor_fiber_time(~con_mice',:),[],1));

subplot(3,4,1)
[~,y1_I] = sort(mean(y1_time(:,time > 0 & time <= 2),2));
imagesc(time,1:size(y1_time,1),y1_time(y1_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'xtickLabel',{'-5','0','2','5','10'},'fontsize',6);
set(gca,'ytick',[1,size(y1_time,1)],'ytickLabel',{'1',size(y1_time,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8); ylabel(('# of Animal-odor pairs'),'fontsize',8);
caxis([-0.2,1]);
    
subplot(3,4,2)
[~,y2_I] = sort(mean(y2_time(:,time > 0 & time <= 2),2));
imagesc(time,1:size(y2_time,1),y2_time(y2_I,:));
set(gca,'xtick',[-4.99 0 2 5 10],'xtickLabel',{'-5','0','2','5','10'},'fontsize',6);
set(gca,'ytick',[1,size(y2_time,1)],'ytickLabel',{'1',size(y2_time,1)},'fontsize',8); 
set(gca,'tickdir','out','ticklength',[0,0],'linewidth',0.0001);
xlabel('Time (s)','fontsize',8);
caxis([-0.2,1]);
    
subplot(3,4,3)
animal_odor_fiber_time = cellfun(@(x) mean(x,1),animal_odor_fiber,'UniformOutput',false);

y1_time = cell2mat(reshape(animal_odor_fiber_time( con_mice',:),[],1));
y2_time = cell2mat(reshape(animal_odor_fiber_time(~con_mice',:),[],1));


hold on;
plot_shadow(time,mean(y1_time,1),sem(y1_time,1),'r');
plot_shadow(time,mean(y2_time,1),sem(y2_time,1),'b');
set(gca,'tickdir','out','ticklength',[0.03 0.03],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlabel('Time (s)','fontsize',10); ylabel(('ΔF/F (%)'),'fontsize',6);
set(gca,'xtick',[-4.99 0 2 5 10],'xtickLabel',{'-5','0','2','5','10'},'fontsize',6);
ylim([-0.2,0.4]);

subplot(3,4,4)
animal_odor_fiber_mean = cell2mat(cellfun(@(x) mean2(x(:,time > 0 & time <= 2)),animal_odor_fiber,'UniformOutput',false));

hold on
y1 = reshape(animal_odor_fiber_mean( con_mice',:),[],1);
y2 = reshape(animal_odor_fiber_mean(~con_mice',:),[],1);
unpaired_stats = unpaired_test({y1,y2});

bar_scatter_alpha(1,y1,'r');
bar_scatter_alpha(2,y2,'b');
title({unpaired_stats.method;[unpaired_stats.stat,', ',unpaired_stats.p]},'FontSize',6);
ylabel(('ΔF/F (%)'),'fontsize',6);
ylim([-0.5 1]);
