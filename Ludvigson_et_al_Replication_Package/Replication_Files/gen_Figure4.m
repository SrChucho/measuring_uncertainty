% -------------------------------------------------------------------------
% Plot Figure 4 in the Paper (Distribution of BYF and BYM)
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
%QQ =[];
%NQ= 1.5*10^(6); %Number of Random Rotation
%randn('state',123456) %Fixed seeds
%parfor iii = 1:NQ
%    %generate a random rotation
%    v = randn(3,3); %First Generate a totally random 3x3 matrix v
%    [q, r]=qr(v,0); %Obtain the orthonormal matrix of v, call it q
%    QQ(:,:,iii)=q*diag(sign(diag(r)));    % this makes the diagonal of r positive
%end
%% Compute IRF

% 5 specifications
% ind_spec == 1: Cov Only -- stored in B_store_hist (No restrictions, all 1.5 million)
% ind_spec == 2: Cov + Event constraints
% ind_spec == 3: Cov + Event + External Variable Constraints
% ind_spec == 4: Cov + External Variable + Event Constraints w/o gbar2
% ind_spec == 5: Cov + Event + External Variable w/o gc1

lo_combine = [];
hi_combine = [];
for ind_spec = 2:5
    
    Ebar = [4.1634,4.0475,4.5672,4.7314];
    
    % Generate B
    B_out= gen_B_hist(vec_A, QQ, X_var, X_exp,const_var,data,dates_lag,lag,Ebar,ind_spec);
    
    if ind_spec == 1
        B_out_cov = B_out;
    elseif ind_spec == 2
        B_out_event = B_out;
    elseif ind_spec == 3
        B_out_both = B_out;
    elseif ind_spec == 4
        B_out_noGFC = B_out;
    elseif ind_spec == 5
        B_out_nostock = B_out;
    end
end

%% Upper Panel

% Load all pre-saved 1.5 million rotations of B
load B_store_hist B_store_cov

B_store_both = reshape(B_out_both,9,size(B_out_both,3))';
B_store_event = reshape(B_out_event,9,size(B_out_event,3))';

scale = 100;
kk_list = [2,8];
plot_list = [1,4,2,5,3,6];
counter = 0;
qc = 0;
fontsize_label = 10;
for kkk = kk_list
    qc = qc+1;
    counter = counter + 1;
    subplot(3,2,qc);
    histogram(B_store_cov(:,kk_list(counter))*scale,15,'Normalization','probability')
    hold on
    histogram(B_store_event(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','r')
    hold on
    histogram(B_store_both(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','k')
    if kkk == 1
        title('$B_{MM}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 2
        title('$B_{YM}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z}=0$, $\bar{g}_{E}\geq 0$', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        
        set(leg,'fontsize',8,'box','off','location','northwest','interpreter','latex')
    elseif kkk == 3
        title('$B_{FM}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 4
        title('$B_{MY}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 5
        title('$B_{YY}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 6
        title('$B_{FY}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 7
        title('$B_{MF}$','fontsize',fontsize_label,'interpreter','latex')
    elseif kkk == 8
        title('$B_{YF}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z}=0$, $\bar{g}_{E}\geq 0$', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        
        set(leg,'fontsize',8,'box','off','location','northeast','interpreter','latex')
    elseif kkk == 9
        title('$B_{FF}$','fontsize',fontsize_label,'interpreter','latex')
    end
    set(gca,'fontsize',fontsize_label)
end


%% Middle Panel

B_store_noGFC = reshape(B_out_noGFC,9,size(B_out_noGFC,3))';


scale = 100;
kk_list = [2,8];
counter = 0;
for kkk = kk_list
    qc = qc+1;
    counter = counter + 1;
    subplot(3,2,qc);
    histogram(B_store_cov(:,kk_list(counter))*scale,15,'Normalization','probability')
    hold on
    histogram(B_store_noGFC(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','r')
    hold on
    histogram(B_store_both(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','k')
    if kkk == 1
        title('$B_{MM}$','fontsize',15,'interpreter','latex')
    elseif kkk == 2
        title('$B_{YM}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z}=0,$ $\bar{g}_{C}\geq 0$, $\bar{g}_{Ej}\geq 0$, $\forall j\neq 2$ ', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        set(leg,'fontsize',8,'box','off','location','northwest','interpreter','latex')
    elseif kkk == 3
        title('$B_{FM}$','fontsize',15,'interpreter','latex')
    elseif kkk == 4
        title('$B_{MY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 5
        title('$B_{YY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 6
        title('$B_{FY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 7
        title('$B_{MF}$','fontsize',15,'interpreter','latex')
    elseif kkk == 8
        title('$B_{YF}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z}=0,$ $\bar{g}_{C}\geq 0$, $\bar{g}_{Ej}\geq 0$, $\forall j\neq 2$ ', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        
        
        set(leg,'fontsize',7.5,'box','off','position',[0.7156,0.5606,0.1869,0.0573],'interpreter','latex')
    elseif kkk == 9
        title('$B_{FF}$','fontsize',15,'interpreter','latex')
    end
    set(gca,'fontsize',fontsize_label)
end


%% Bottom Panel

B_store_nostock = reshape(B_out_nostock,9,size(B_out_nostock,3))';

scale = 100;
kk_list = [2,8];
counter = 0;
for kkk = kk_list
    qc = qc+1;
    counter = counter + 1;
    subplot(3,2,qc);
    histogram(B_store_cov(:,kk_list(counter))*scale,15,'Normalization','probability')
    hold on
    histogram(B_store_nostock(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','r')
    hold on
    histogram(B_store_both(:,kk_list(counter))*scale,15,'Normalization','probability','FaceColor','k')
    
    if kkk == 1
        title('$B_{MM}$','fontsize',15,'interpreter','latex')
    elseif kkk == 2
        title('$B_{YM}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z} = 0$, $\bar{g}_{Cj}\geq 0$, $\forall j\neq 1$, $\bar{g}_{E} \geq 0$ ', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        set(leg,'fontsize',8,'box','off','location','northwest','interpreter','latex')
        
        
        set(leg,'fontsize',8,'box','off','location','northwest','interpreter','latex')
    elseif kkk == 3
        title('$B_{FM}$','fontsize',15,'interpreter','latex')
    elseif kkk == 4
        title('$B_{MY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 5
        title('$B_{YY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 6
        title('$B_{FY}$','fontsize',15,'interpreter','latex')
    elseif kkk == 7
        title('$B_{MF}$','fontsize',15,'interpreter','latex')
    elseif kkk == 8
        title('$B_{YF}$','fontsize',fontsize_label,'interpreter','latex')
        leg = legend('$\bar{g}_{Z} = 0$','$\bar{g}_{Z} = 0$, $\bar{g}_{Cj}\geq 0$, $\forall j\neq 1$, $\bar{g}_{E} \geq 0$ ', ...
            '$\bar{g}_{Z}=0$, $\bar{g}_{C}\geq 0,$ $\bar{g}_{E}\geq 0$');
        set(leg,'fontsize',7.5,'box','off','location','northeast','interpreter','latex')
    elseif kkk == 9
        title('$B_{FF}$','fontsize',15,'interpreter','latex')
    end
    set(gca,'fontsize',fontsize_label)
end

set(gcf,'paperpositionmode','manual','paperunits','inches');
dim = [10,8];
set(gcf,'papersize',dim,'paperposition',[0,0,dim]);
print(gcf,'-dpdf','./Figure/Figure4.pdf');