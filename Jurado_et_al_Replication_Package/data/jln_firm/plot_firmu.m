% -------------------------------------------------------------------------
% Plot aggregate uncertainty estimates
% -------------------------------------------------------------------------

% Load estimates
clear; clc; close all;
load ferrors;
load firmu;
T = length(dates);

% Load raw industrial production
[data,txt] = xlsread('jlnrawdata.xlsx',1);
ipg = data(3:end,6);
for i = 1:3:length(ipg)-1;
    qipg(1+(i-1)/3,:) = mean(ipg(i:i+2,:));
end;
ipg = qipg;
ipg = [NaN;log(ipg(2:end)./ipg(1:end-1))];
ipg = tsmovavg(ipg,'s',4,1);
ipg = ipg(end-T+1:end)*100;

% Plot csa estimates
fig = figure(1);
gr  = [0,0.5,0];
u1  = utcsa(:,1);
u2  = utcsa(:,2);
u4  = utcsa(:,4);
pt  = std(zscore(xt),0,2);
pt  = pt(end-T+1:end);
plot(dates,pt,'color',gr,'linewidth',0.5); hold on;
plot(dates,u1,'b','linewidth',1.5);
plot(dates,u2,'--k','linewidth',1.5);
plot(dates,u4,'-.r','linewidth',1.5);
plot([dates(1),dates(end)],(mean(u1)+1.65*std(u1)).*[1,1],'--b');
plot([dates(1),dates(end)],(mean(u2)+1.65*std(u2)).*[1,1],'--k');
plot([dates(1),dates(end)],(mean(u4)+1.65*std(u4)).*[1,1],'--r');
txt = '$\overline{\mathcal{U}}^y_t(4)$';
text(1976.5,1.1,txt,'interpreter','latex','color','r');
txt  = '$\overline{\mathcal{U}}^y_t(2)$';
text(1984,0.33,txt,'interpreter','latex','color','k');
txt  = '$\overline{\mathcal{U}}^y_t(1)$';
text(2004,0.54,txt,'interpreter','latex','color','b');
txt  = '$\mathcal{D}^B_t$';
text(1996,1.6,txt,'interpreter','latex','color',gr);
lab1 = sprintf('h = 1, corr with IP = %0.2f',corr(u1,ipg));
lab2 = sprintf('h = 2, corr with IP = %0.2f',corr(u2,ipg));
lab3 = sprintf('h = 4, corr with IP = %0.2f',corr(u4,ipg));
lab4 = sprintf('Dispersion firm profits, corr with IP = %0.2f',corr(pt,ipg));
annotation('arrow',[0.37,0.28],[0.17,0.25],'linewidth',0.5);
annotation('arrow',[0.28,0.28],[0.41,0.31],'linewidth',0.5,'color','r');
leg = legend(lab4,lab1,lab2,lab3);
set(leg,'location','northwest','box','off');
xlim([dates(1),dates(end)]);
ylim([0.4,1.4]);
ylim([0.2,2.4]);
rshade(dates);

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','firmu_csa');

% Plot pca estimates
fig = figure(2);
u1  = utpca(:,1);
u2  = utpca(:,2);
u4  = utpca(:,4);
plot(dates,u1,'b','linewidth',1.5); hold on;
plot(dates,u2,'k','linewidth',1);
plot(dates,u4,'-.r','linewidth',1.5);
plot([dates(1),dates(end)],(mean(u1)+1.65*std(u1)).*[1,1],'--b');
plot([dates(1),dates(end)],(mean(u2)+1.65*std(u2)).*[1,1],'--k');
plot([dates(1),dates(end)],(mean(u4)+1.65*std(u4)).*[1,1],'--r');
txt4 = '$\widehat{\mathcal{U}}^y_t(4)$';
txt2  = '$\widehat{\mathcal{U}}^y_t(2)$';
txt1  = '$\widehat{\mathcal{U}}^y_t(1)$';
text(1977.5,0.81,txt4,'interpreter','latex','color','r');
text(1983,0.5,txt2,'interpreter','latex','color','k');
text(2004,0.55,txt1,'interpreter','latex','color','b');
lab1 = sprintf('h = 1, corr with IP = %0.2f',corr(u1,ipg));
lab2 = sprintf('h = 2, corr with IP = %0.2f',corr(u2,ipg));
lab3 = sprintf('h = 4, corr with IP = %0.2f',corr(u4,ipg));
annotation('arrow',[0.36,0.3],[0.2,0.25],'linewidth',0.5);
leg = legend(lab1,lab2,lab3);
set(leg,'location','northwest','box','off');
xlim([dates(1),dates(end)]);
ylim([0.4,1.4]);
rshade(dates);

% Print figure
dim = [6,5];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(fig,'-dpdf','firmu_pca');