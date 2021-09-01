cd('E:\Documents\MATLAB\2020-07-02 OB-Syn-MT-Array Fengjiao Chen\02_Mat\SPKC_Example');

load('con_07_SPKC02.mat');
load('syn_09_SPKC03.mat');

Ctrl = SPKC02;
Syn  = SPKC03;

fs = 1/SPKC03_ts_step;

%% Plot
time = 10:14;

fig1 = figure('name','Spikes_Raw_Trace','numbertitle','off');
set(gcf,'unit','centimeters','PaperUnits','centimeters','PaperPosition',[0,0,20,15],'color','w','PaperSize',[20,15]);
subplot(3,2,1);
plot( Ctrl(time(1)*fs:time(end)*fs),'color','k','LineWidth',0.5);
box off
set(gca,'tickdir','out','ticklength',[0.03 0.03],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
ylim([-0.2 0.2]);

subplot(3,2,2);
plot( Syn(time(1)*fs:time(end)*fs),'color','r','LineWidth',0.5);
box off
set(gca,'tickdir','out','ticklength',[0.03 0.03],...
    'linewidth',0.75,'Fontname', 'Arial','fontsize',6);
ylim([-0.4 0.4]);