% Set default parameters
global fs bin pre post stimu
fs = 1000; bin = 0.1;
pre = 4;   post = 6;
stimu = 2;

% Select the folders to be analyzed
FolderList = uipickfiles;

%% Pre-analysis
unit_odor_mfr = [];
unit_mice_type = [];
for i = 1:size(FolderList,2)
    
    % Determine the current subfolder name
    selected_path = FolderList{i};
    cd(selected_path);
    work_folder = selected_path(find(selected_path == '\',1,'last')+1:end);
    disp([work_folder,' is being processed . . .']);
    
    % Load 'Odor¡¢Behavior' raw data from 'Behavior.txt'
    change_suffix('dat','txt');
    try
        odor = str2num(cell2mat(load_result('Behavior.txt',4)));
    catch
        odor = textread('Behavior.txt','%n');
    end
    
    % Find event
    event = find_event(2000);

    % Load Spikes
    Spikes = load('Spikes.mat'); Spikes = rmfield(Spikes,'spikes_all');
    
    % Check number of Trials and Behaviors
    if length(event) ~= length(odor)
        error('Error>> The number of Trials is different from the number of Behaviors !')
    end
    
    % Calculate the Mean firing rate in each Trial £¨trial# - time - unit#£©
    spikes_mfr = structfun(@(x) firing_rate(x,event)./bin,Spikes,'UniformOutput',false);
    
    % Sort Period | Odor
    % odor_period data structure : Odor * Period
    unit_names = fieldnames(spikes_mfr);
    for unit_index = 1:size(unit_names,1)
        unit_mfr = eval(['spikes_mfr.',unit_names{unit_index},';']);
        for odor_chan = min(odor):max(odor)
            odor_mfr{:,odor_chan+1} = unit_mfr(odor == odor_chan,:);
        end
        unit_odor_mfr  = [unit_odor_mfr; odor_mfr];
        unit_mice_type = [unit_mice_type;contains(work_folder,'Ctrl')]; % 1 --> CON; 0 --> EXP
    end

    clearvars -except fs pre post bin stimu FolderList i unit_odor_mfr unit_mice_type
end
unit_mice_type = repmat(unit_mice_type,1,8);

%% Analysis
% Load mat file
load('E:\Documents\MATLAB\2020-07-02 OB-Syn-MT-Array Fengjiao Chen\02_Mat\Spikes_Features_Synuclein.mat');

prebase    = 0.2;
basetime   = (pre-stimu-prebase)/bin+1:(pre-prebase)/bin;
evokedtime = (pre/bin +1:(pre+ stimu)/bin);

% Respones
evoked_response = cell2mat(cellfun(@(x) spikes_response(x),unit_odor_mfr(:,:),'UniformOutput',false));

% Spotaneous MFR
spotaneous_mfr = mean(cell2mat(cellfun(@(x) mean2(x(:,basetime)),unit_odor_mfr,'UniformOutput',false)),2);

% Evoked  MFR
evoked_mfr = cell2mat(cellfun(@(x) mean2(x(:,evokedtime)),unit_odor_mfr,'UniformOutput',false))-spotaneous_mfr;

% SNR
evoked_snr = evoked_mfr./spotaneous_mfr;

%% Plot spotaneous_mfr
fig1 = figure('name','OB-Syn_SpontaneousEvoked_Spikes','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);

y1 = spotaneous_mfr(unit_mice_type(:,1)==1);
y2 = spotaneous_mfr(unit_mice_type(:,1)==0);

% Histogram of spotaneous_mfr
subplot(3,3,1)
hold on
BinWidth = 2;
histogram(y1,'BinWidth',BinWidth,'LineWidth',0.005,'FaceColor','r','FaceAlpha',0.5);
histogram(y2,'BinWidth',BinWidth,'LineWidth',0.005,'FaceColor','b','FaceAlpha',0.5);
box off
set(gca,'tickdir','out','ticklength',[0.03 0.03],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlabel('Mean firing rate (Hz)','fontsize',8);
ylabel('Unit count','fontsize',8);
xlim([0 50]);

% Cumulative fraction units of spotaneous_mfr
[~,p,ks2stat] = kstest2(y1,y2);

subplot(3,3,2)
hold on
h1 = cdfplot(y1);
h2 = cdfplot(y2);
set( h1, 'LineStyle', '-', 'Color', 'r', 'linewidth', 0.75);
set( h2, 'LineStyle', '-', 'Color', 'b', 'linewidth', 0.75);
box off; title(''); grid off;
set(gca,'tickdir','out','ticklength',[0.025 0.025],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlim([0 50]);
xlabel('Mean firing rate (Hz)','fontsize',8);
ylabel('Cumulative fraction units','fontsize',8);
title({'Two-sample K-S test';[['Z', ' = ',num2str(ks2stat)],', ',['P = ',num2str(p)]]},'FontSize',6);

% Bar_scatter of spotaneous_mfr
unpaired_stats = unpaired_test({y1,y2});

subplot(3,3,3)
hold on
bar_scatter_alpha(1,y1,'r');
bar_scatter_alpha(2,y2,'b');
ylabel('Mean firing rate (Hz)','fontsize',8);
ylim([0 50])
title({unpaired_stats.method;[unpaired_stats.stat,', ',unpaired_stats.p]},'FontSize',6);


%% Plot evoked_mfr
con_response = [numel(evoked_mfr(unit_mice_type==1 & evoked_response==1));...
    numel(evoked_mfr(unit_mice_type==1 & evoked_response==2));...
    numel(evoked_mfr(unit_mice_type==1 & evoked_response==0))];
exp_response = [numel(evoked_mfr(unit_mice_type==0 & evoked_response==1));...
    numel(evoked_mfr(unit_mice_type==0 & evoked_response==2));...
    numel(evoked_mfr(unit_mice_type==0 & evoked_response==0))];

% stackedbar of unit response
subplot(3,3,4)
h = bar([con_response/sum(con_response), exp_response/sum(exp_response)]'*100, 0.4, 'stacked');
set(h,{'FaceColor'},{color('orange'); color('green'); color('grey')});
box off
set(gca,'tickdir','out','ticklength',[0.03 0.03],'xtick',[],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
ylabel('Percent of response','fontsize',8);

% Bar_scatter of evoked_mfr
con_exci_mfr = evoked_mfr(unit_mice_type==1 & evoked_response==1);
con_inhi_mfr = evoked_mfr(unit_mice_type==1 & evoked_response==2);
exp_exci_mfr = evoked_mfr(unit_mice_type==0 & evoked_response==1);
exp_inhi_mfr = evoked_mfr(unit_mice_type==0 & evoked_response==2);

subplot(3,3,5)
hold on
y1 = con_exci_mfr;
y2 = exp_exci_mfr;

unpaired_stats = unpaired_test({y1,y2});

bar_scatter_alpha(1,y1,'r');
bar_scatter_alpha(2,y2,'b');
ylabel('Mean firing rate (Hz)','fontsize',8);
ylim([0 30]);
title({unpaired_stats.method;[unpaired_stats.stat,', ',unpaired_stats.p]},'FontSize',6);

subplot(3,3,6)
hold on
y1 = con_inhi_mfr;
y2 = exp_inhi_mfr;

unpaired_stats = unpaired_test({y1,y2});

bar_scatter_alpha(1,y1,'r');
bar_scatter_alpha(2,y2,'b');
ylabel('Mean firing rate (Hz)','fontsize',8);
ylim([-15 0])
title({unpaired_stats.method;[unpaired_stats.stat,', ',unpaired_stats.p]},'FontSize',6);

%% Plot evoked_snr
%{ 
% Histogram of spotaneous_mfr
con_snr = evoked_snr(unit_mice_type==1 & evoked_response~=0);
exp_snr = evoked_snr(unit_mice_type==0 & evoked_response~=0);

y1 = con_snr;
y2 = exp_snr;
subplot(3,3,7)
hold on
BinWidth = 0.05;
histogram(y1,'BinWidth',BinWidth,'LineWidth',0.005,'FaceColor','r','FaceAlpha',0.5);
histogram(y2,'BinWidth',BinWidth,'LineWidth',0.005,'FaceColor','b','FaceAlpha',0.5);
box off
set(gca,'tickdir','out','ticklength',[0.03 0.03],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlabel('SNR','fontsize',8);
ylabel('Unit count','fontsize',8);
ylim([0 60])

% Cumulative fraction units of spotaneous_mfr
subplot(3,3,8)
[~,p,ks2stat] = kstest2(y1,y2);
hold on
h1 = cdfplot(y1);
h2 = cdfplot(y2);
set( h1, 'LineStyle', '-', 'Color', 'r', 'linewidth', 1);
set( h2, 'LineStyle', '-', 'Color', 'b', 'linewidth', 1);
ylabel('Cumulative fraction units','fontsize',8);
box off; title(''); grid off;
set(gca,'tickdir','out','ticklength',[0.025 0.025],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
xlabel('SNR','fontsize',8);
z_value = ['\itz', ' = ',num2str(ks2stat)];
p_value = ['\itp = ',num2str(p)];
title([z_value,',   ',p_value],'FontSize',8);
%}

%% Errorbar comparing evoked-MFR by different odors
%{
subplot(3,3,9)
y1 = evoked_mfr(unit_mice_type(:,1)==1,:);
y2 = evoked_mfr(unit_mice_type(:,1)==0,:);

% ANOVAN
y = reshape([y1;y2],[],1);
g1_0 = [ones(size(y1)); zeros(size(y2))];
for odor_chan = 1:8
    g2_0(:,odor_chan) = ones(size(g1_0,1),1)*odor_chan;
end
g1 = reshape(g1_0,[],1);
g2 = reshape(g2_0,[],1);
[~,table,stats,~]= anovan(y,{g1,g2},'display','off');
 unpaired_stats.method = 'Two-way ANOVA';

unpaired_stats.table = table;
unpaired_stats.stat1 = [table{1,6},'(',num2str(table{2,3}),',',num2str(table{4,3}),')', ' = ', num2str(table{2,6})];
unpaired_stats.stat2 = [table{1,6},'(',num2str(table{3,3}),',',num2str(table{4,3}),')', ' = ', num2str(table{3,6})];
unpaired_stats.p1 = ['P = ',num2str(table{2,7})];
unpaired_stats.p2 = ['P = ',num2str(table{3,7})];

c_val = multcompare(stats,'Display','off');
row = find(c_val(:,6)<0.05);
for r = 1:length(row)
    unpaired_stats.posthoc{r,:} = [num2str(c_val(row(r),1)),' vs ' ,num2str(c_val(row(r),2)),', P = ',num2str(c_val(row(r),6))];
end

hold on
errorbar(1:8,mean(y1,1),sem(y1,1),'-or','MarkerSize',2,'MarkerFacecolor','r','CapSize',3)
errorbar(1:8,mean(y2,1),sem(y2,1),'-ob','MarkerSize',2,'MarkerFacecolor','b','CapSize',3)

box off
set(gca,'tickdir','out','ticklength',[0.02 0.02],...
    'linewidth',0.5,'Fontname', 'Arial','fontsize',6);
xlabel('Odortype','fontsize',8);
ylabel('¦¤ Mean firing rate (Hz)','fontsize',8);
xlim([0.5 8.5]);
ylim([-20 20]);
title({unpaired_stats.method;[unpaired_stats.stat1,', ',unpaired_stats.p1];...
    [unpaired_stats.stat2,', ',unpaired_stats.p2]},'FontSize',6);
disp(unpaired_stats.posthoc);

%}
