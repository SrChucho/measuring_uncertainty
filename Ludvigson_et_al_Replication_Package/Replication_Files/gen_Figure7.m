% -------------------------------------------------------------------------
% Plot Figure 7 in the Paper (Economic Policy Uncertainty)
% -------------------------------------------------------------------------

restoredefaultpath
clear; clc; close all;
addpath './supporting_code'

%% Load Data
IRF_lim = 60; % IRF length
lag = 6;      % VAR Lag

% Load Data
U_data = importdata('Replication_Data.xlsx');
U_data = U_data.data.Data;
dates = 1960 + 7/12 : 1/12 : 2015+4/12;
Uf = U_data(:,1);
Um = U_data(:,2);
ip = U_data(:,4);
Vs = U_data(:,6);
Gold_o = U_data(:,7);
Gold_o = [NaN;diff(log(Gold_o))];

% Load EPU
EPU = U_data(:,8);
EPU_news =  U_data(:,9);

%% Trim Sample (Use EPU sample)
ind = isnan(EPU) ==0;
Um = Um(ind);
Uf = Uf(ind);
ip = ip(ind);
Vs = Vs(ind);
Gold_o = Gold_o(ind);
dates = dates(ind);
dates_lag = dates(lag+1:end);

% Standardize EPU
EPU = zscore(EPU(ind))*std(Um) + mean(Um);
EPU_news = zscore(EPU_news(ind))*std(Um) + mean(Um);


% 2 Specifications
% k == 1: EPU Base
% k == 2: EPU News

for k = 1:2
    %% Stack variables
    if k == 1
        X = [EPU,ip,Uf];
    else
        X = [EPU_news,ip,Uf];
    end
    X_var = X(lag+1:end,:);
    X_exp = [];
    for jjj = 1:lag
        X_exp= [X_exp, X(lag+1-jjj:end-jjj,:)];
    end
    
    %% Reduced Form VAR
    T = size(X_var,1);
    data = [Vs(lag+1:end), Gold_o(lag+1:end)];
    reg31 = regstats(X_var(:,1), X_exp);
    reg32 = regstats(X_var(:,2), X_exp);
    reg33 = regstats(X_var(:,3), X_exp);
    eta_m_OLS = reg31.r;
    eta_y_OLS = reg32.r;
    eta_f_OLS = reg33.r;
    eta = [eta_m_OLS,eta_y_OLS,eta_f_OLS];
    eta_vec_OLS = eta;
    vec_A = [reg31.beta(2:end);reg32.beta(2:end);reg33.beta(2:end)];
    const_var = [reg31.beta(1), reg32.beta(1), reg33.beta(1)];
    
    %% Generate Random Rotation
    load QQ_full % Load pre-saved seeds used in the paper
    
    %Comment out to generate new random rotations
    %QQ =[];
    %NQ= 1.5*10^(6); %Number of Random Rotation
    %parfor iii = 1:NQ
    %    %generate a random rotation
    %    v = randn(3,3); %First Generate a totally random 3x3 matrix v
    %    [q, r]=qr(v,0); %Obtain the orthonormal matrix of v, call it q
    %    QQ(:,:,iii)=q*diag(sign(diag(r)));    % this makes the diagonal of r positive
    %end
    %% Compute e
    Ebar = [5.4261,3.8529,2];
    % Generate B
    B_out= gen_B_EPU(vec_A, QQ, X_var, X_exp,const_var,data,dates_lag,lag,Ebar);
    % Gen IRF
    [lo,hi] =gen_lo_hi(vec_A, B_out,X_var, X_exp,const_var,data,lag);
    lo_combine(:,:,k) = lo(:,:);
    hi_combine(:,:,k) = hi(:,:);
    
end
%% Plot Time Series (Left Panel)

fig = figure;
f1 = plot(dates, zscore(EPU(:,1)), 'color','r', 'Linewidth',2);
hold on 
f2 = plot(dates, zscore(EPU_news(:,1)), 'color','b', 'Linewidth',1.5);
hold on 
plot(dates, 1.65*ones(length(dates),1), 'k--', 'Linewidth',1)
xlabel('Year','interpreter','latex','Fontsize',15)
ylabel('Standardized unit','fontsize',15)
title('Policy Uncertainty','interpreter','latex','Fontsize',15)
xlim([dates(1),dates(end)]);
rshade(dates);
leg=legend([f1,f2],'EPU','EPN');
set(leg,'fontsize',15,'location','northwest','box','off')

dim = [8,8];
set(gcf,'paperpositionmode','manual','paperunits','inches');
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure7_left');


%% Plot IRF (Right Panel)
ylab{1}=['EPU'];
ylab{2}=['Y'];
ylab{3}=['U_{F}'];
fig= plot_irf_fan_two_EPU(hi_combine(:,:,[1,2,2]),lo_combine(:,:,[1,2,2]),IRF_lim,ylab);
set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [8,6];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure7_right');
