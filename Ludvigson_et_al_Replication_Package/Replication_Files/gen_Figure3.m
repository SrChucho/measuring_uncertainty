% -------------------------------------------------------------------------
% Plot Figure 3 in the Paper (IRFs under Di¤erent Parameterizations of Big Shocks)
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

%% Trim Sample (If needed)
ind = isnan(Um) ==0;
Um = Um(ind);
Uf = Uf(ind);
ip = ip(ind);
Vs = Vs(ind);
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
load QQ_full
% QQ =[];
% NQ= 1.5*10^(6); %Number of Random Rotation
% randn('state',123456) %Fixed seeds
% parfor iii = 1:NQ
%     %generate a random rotation
%     v = randn(3,3); %First Generate a totally random 3x3 matrix v
%     [q, r]=qr(v,0); %Obtain the orthonormal matrix of v, call it q
%     QQ(:,:,iii)=q*diag(sign(diag(r)));    % this makes the diagonal of r positive
% end
%% Compute IRF

% 5 specifications
% ind_spec == 1: 1987 and 1970 at 75th percentile, others fixed at 75th
% ind_spec == 2: 1987 and 1970 at 85th percentile, others fixed at 75th
% ind_spec == 3: 1987 and 1970 at median, others fixed at 75th
% ind_spec == 4: Lehman event  at 85th percentile, others fixed at 75th
% ind_spec == 5: Lehman event  at median , others fixed at 75th

lo_combine = [];
hi_combine = [];
for ind_spec = 1:5
    if ind_spec == 1
       Ebar = [4.1634,4.0475,4.5672,4.7314];
    elseif ind_spec == 2
       Ebar = [4.4415,4.3915,4.5672,4.7314]; 
    elseif ind_spec == 3
       Ebar = [2.7744,2.3317,4.5672,4.7314]; 
    elseif ind_spec == 4
       Ebar = [4.1634,4.0475,4.9305,5.0442];
    elseif ind_spec == 5
       Ebar = [4.1634,4.0475,1.9461,2.3317];
    end
    
    % Generate B
    B_out= gen_B(vec_A, QQ, X_var, X_exp,const_var,data,dates_lag,lag,Ebar);
    % Gen IRF
    [lo,hi]=gen_lo_hi(vec_A, B_out,X_var, X_exp,const_var,data,lag);
    
    lo_combine(:,:,ind_spec ) = lo;
    hi_combine(:,:,ind_spec ) = hi;
end

%% Plot

% Left Panel
for k = 1:3
lo_lef(:,:,k) = lo_combine(:,:,k);
hi_lef(:,:,k) = hi_combine(:,:,k);
end
ylab{1}=['U_{M}'];
ylab{2}=['Y'];
ylab{3}=['U_{F}'];
fig= plot_irf_fan_three(hi_lef,lo_lef,IRF_lim,ylab);
set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [10,8];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure3_left');

% Right Panel
counter = 0;
for k = [1,4,5]
    counter = counter +1;
lo_right(:,:, counter ) = lo_combine(:,:,k);
hi_right(:,:, counter ) = hi_combine(:,:,k);
end
ylab{1}=['U_{M}'];
ylab{2}=['Y'];
ylab{3}=['U_{F}'];
fig= plot_irf_fan_three(hi_right,lo_right,IRF_lim,ylab);
set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [10,8];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure3_right');
