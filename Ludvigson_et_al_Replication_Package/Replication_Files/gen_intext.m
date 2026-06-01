% -------------------------------------------------------------------------
% Generate In-text numbers
% -------------------------------------------------------------------------

restoredefaultpath
clear; close all;
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
CRSP = U_data(:,6);
Gold_o = U_data(:,7);
Gold_o = [NaN;diff(log(Gold_o))];
VIX = U_data(:,12);
%% Trim Sample (If needed)
ind = isnan(Um) ==0;
Um = Um(ind);
Uf = Uf(ind);
ip = ip(ind);
CRSP = CRSP(ind);
Gold_o = Gold_o(ind);
dates = dates(ind);
dates_lag = dates(lag+1:end);

%% Stack variables
X = [Um,ip,Uf];
X_var = X(lag+1:end,:);
X_exp = [];
for jjj = 1:lag
    X_exp= [X_exp, X(lag+1-jjj:end-jjj,:)];
end

%% Reduced Form VAR
T = size(X_var,1);
reg31 = regstats(X_var(:,1), X_exp);
reg32 = regstats(X_var(:,2), X_exp);
reg33 = regstats(X_var(:,3), X_exp);
eta_m_OLS = reg31.r;
eta_y_OLS = reg32.r;
eta_f_OLS = reg33.r;
eta = [eta_m_OLS,eta_y_OLS,eta_f_OLS];
vec_A = [reg31.beta(2:end);reg32.beta(2:end);reg33.beta(2:end)];
const_var = [reg31.beta(1), reg32.beta(1), reg33.beta(1)];

%% Generate Random Rotation
load QQ_full % Load pre-saved seeds used in the paper

dates0 = dates(1+lag:end); %Trim the sample based on VAR lag
OMEGA= cov(eta);
BB0=chol(OMEGA,'lower'); %Generate lower cholesy factorization of covariance matrix eta

% Generate Dates
index_87 = [find(abs(dates0 - (1987 + 10/12))<10^(-6))]; %Find 1987:10
index_70 = [find(abs(dates0 - (1970 + 12/12))<10^(-6))]; %Find 1970:12
index_089 = [find(abs(dates0 - (2008 + 9/12))<10^(-6))]; %Find 2008:09

% Generate e
parfor iii = 1:size(QQ,3)  
    % rotate and normalize B
    B_initial=BB0*QQ(:,:,iii);
    B_est=B_initial*diag(sign(diag(B_initial))); %normalize so that self-irfs are positive
    
    ehat = eta(:,:)*inv(B_est)';
    ev_store(iii,:) = [ehat(index_87,3),ehat(index_089,3),ehat(index_70,1),ehat(index_089,1),];
end

%% Page 10, Line 12

VIX_sample = logical(1-isnan(VIX)); %VIX sample starts from 1990:01
disp('Page 10, Line 12')
disp('Correlation between VIX and Log Returns')
disp(corr(VIX(VIX_sample), log(1+(CRSP(VIX_sample))/100)))

%% Page 22, Lines 11-14

disp('Page 22, Lines 11-14')
disp('75th percentile of eF at 1987:10')
disp(prctile(ev_store(:,1),75))
disp('75th percentile of eF at 2008:09')
disp(prctile(ev_store(:,2),75))
disp('75th percentile of eM at 1970:12')
disp(prctile(ev_store(:,3),75))
disp('75th percentile of eM at 2008:09')
disp(prctile(ev_store(:,4),75))
disp('If shocks were Gaussian, the probability of a shock of 4 std is ')
disp(normpdf(4))


