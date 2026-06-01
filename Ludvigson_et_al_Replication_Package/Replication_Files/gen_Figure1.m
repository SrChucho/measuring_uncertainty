% -------------------------------------------------------------------------
% Plot aggregate uncertainty Um and Uf -- Figure 1
% -------------------------------------------------------------------------

restoredefaultpath
addpath './supporting_code'
clear; clc; close all;

%% Load Data
Data_set = importdata('Replication_Data.xlsx'); %Raw Data from 1960:7 to 2015:4
Data_set = Data_set.data.Data; 

% Load Fiancial Uncertainty
utfcsa  = Data_set(:,1);

% Load Macro Uncertainty
utycsa = Data_set(:,2);

% Load Real Uncertainty
U_R = Data_set(:,3);
dates =1960+7/12:1/12 :2015+4/12;

T = length(dates);

% Load ipg
ipg = Data_set(:,5);

%% Plot Time Seris Um and Uf and UR -- Combined
ind = zscore(utycsa(:,1)) > 1.65;
ind_f = zscore(utfcsa(:,1)) > 1.65;
ind_r = zscore(U_R(:,1)) > 1.65;
stan_um = zscore(utycsa(:,1));
stan_uf = zscore(utfcsa(:,1));
stan_ur = zscore(U_R(:,1));

fig = figure;
subplot(2,1,1)
f1=plot(dates, zscore(utycsa(:,1)), 'b', 'Linewidth',1.5);
hold on 
plot(dates(ind), stan_um(ind), 'k*', 'Linewidth',1.5)
hold on
plot(dates, 1.65*ones(length(dates),1), 'k--', 'Linewidth',1)
xlabel('Year','interpreter','latex','Fontsize',15)
text(2009.8,1.9,'1.65 std','color','r');
title('Aggregate Macro Uncertainty $U_{M}$ ','interpreter','latex','Fontsize',15)
lab1 = sprintf('$U_{M}$, corr with IP = %0.2f',corr(utycsa(:,1),ipg))   ;
rshade(dates);
leg = legend(f1,lab1);
set(leg,'location','northwest','box','off','interpreter','latex','fontsize',12);
xlim([dates(1),dates(end)]);

subplot(2,1,2)
f1=plot(dates, zscore(utfcsa(:,1)), 'r', 'Linewidth',1.5);
hold on 
plot(dates(ind_f), stan_uf(ind_f), 'k*', 'Linewidth',1.5)
hold on
plot(dates, 1.65*ones(length(dates),1), 'k--', 'Linewidth',1)
xlabel('Year','interpreter','latex','Fontsize',15)
text(2009.8,1.9,'1.65 std','color','r');
title('Aggregate Financial Uncertainty $U_{F}$ ','interpreter','latex','Fontsize',15)
lab1 = sprintf('$U_{F}$, corr with IP = %0.2f',corr(utfcsa(:,1),ipg))   ;
rshade(dates);
leg = legend(f1,lab1);
set(leg,'location','northwest','box','off','interpreter','latex','fontsize',12);
xlim([dates(1),dates(end)]);
dim = [8,8];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure1');

