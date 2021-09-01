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

     
    % FFT
    window = hamming(pre*fs);
    noverlap = 0;
    nfft = max(256,2^nextpow2(length(window)));
    for k = 1:length(event)
        x = lfp(k,:);
        [~,freq,time,power(:,k)] = spectrogram(x(1,1:pre*fs),window,noverlap,nfft,fs);
    end
    
    % Power to dB
    power = 10*log10(power*10^6);
    
    % Mean Trial
    power_mice(:,f) = mean(power,2);
        
    clearvars -except fs pre post stimu FolderList f freq unit_mice_type power_mice

end

%% Plot
[theta,beta,loga,higa] = find_band(freq);
band{1} = beta; 
band{2} = higa;

y1{1} = power_mice(beta,  unit_mice_type);
y2{1} = power_mice(beta, ~unit_mice_type);
y1{2} = power_mice(higa,  unit_mice_type);
y2{2} = power_mice(higa, ~unit_mice_type);  

z1{1} = mean(power_mice(beta,  unit_mice_type),1);
z2{1} = mean(power_mice(beta, ~unit_mice_type),1);
z1{2} = mean(power_mice(higa,  unit_mice_type),1);
z2{2} = mean(power_mice(higa, ~unit_mice_type),1); 

fig1 = figure('name','LFP_Spontaneous_Synuclein','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);

for b = 1:2
    subplot(3,4,b)
    hold on
    plot_shadow(freq(band{b}),mean(y1{b},2),sem(y1{b},2),'r');
    plot_shadow(freq(band{b}),mean(y2{b},2),sem(y2{b},2),'b');
    box off
    set(gca,'tickdir','out','ticklength',[0.03 0.03],...
        'linewidth',0.75,'Fontname', 'Arial','fontsize',8);
    xlabel('Frequency (Hz)','fontsize',10);ylabel('Power (dB)','fontsize',10);
end

for b = 1:2
    subplot(3,4,b+2)
    unpaired_stats = unpaired_test({z1{b}',z2{b}'});
    hold on
    bar_scatter_alpha(1,z1{b},'r');
    bar_scatter_alpha(2,z2{b},'b');
    ylabel('Power (dB)','fontsize',10);
    title({unpaired_stats.method;[unpaired_stats.stat,', ',unpaired_stats.p]},'FontSize',6);
end






