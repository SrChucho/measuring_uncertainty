function [lo,hi,out_VD,B_allsol,c_out,ehat_maxG] = gen_lo_hi(vec_A, B_in,X_var, X_exp,const_var,data,lag)

% External Variables
Vs = data(:,1);
Gold = data(:,2);

% Parameter
Frac_sol = 100; %Retain 100 percent of solutions
IRF_lim = 60;
NQ = size(B_in,3);

Q = [[vec_A(1:3*lag)';vec_A(3*lag+1:6*lag)';vec_A(6*lag+1:9*lag)'];[eye(3*(lag-1)),zeros(3*(lag-1),3)]];
J = [eye(3),zeros(3,3*(lag-1))];
A_minus_ind = ones(3,1);
eta = (X_var' - repmat(const_var,size(X_var,1),1)'-[vec_A(1:3*lag)';vec_A(3*lag+1:6*lag)';vec_A(6*lag + 1:9*lag)']*X_exp')';


% Solve B
parfor iii = 1:NQ
    
    % rotate and normalize B
    B_est= B_in(:,:,iii);
    ehat = eta(:,:)*inv(B_est)';
    
    % Compute Correlations
    phi1 = corr(Vs, ehat(:,1));
    phi2 = corr(Vs, ehat(:,3));
    phi1_gold = corr(Gold, ehat(:,1));
    phi2_gold = corr(Gold, ehat(:,3));
    
    C_est = sqrt(phi2.^2 + phi1_gold.^2);
    
    % Winnowing constraint
    ind_ef(iii) = 1;
    
    B_store(:,:,iii) = B_est;
    c_store(iii,:) = [phi1,phi2,phi1_gold,phi2_gold,C_est];
    
    G_temp = [0 - phi1,0 - phi2, phi1_gold- 0, phi2_gold- 0];
    G_sr(iii) = sqrt(G_temp*G_temp');
    
end

B_allsol = B_store(:,:,logical(ind_ef));
c_out = c_store(logical(ind_ef),:);
Ns = size(B_allsol,3);
ctr = size(B_allsol,3);

% Max G solution
% It returns the max-G result that satisfy winnowing contraint
[~,max_Ck_indx] = max(G_sr(logical(ind_ef)));
ehat_maxG = eta*inv(B_allsol(:,:,max_Ck_indx))';

for kk = 1:Ns
    
    
    B_0 = B_allsol(:,:,kk);
    % Get the true IRF from data
    for ttt = 0:IRF_lim
        PSI_IRS_m(:,ttt+1) = (J*(Q^ttt)*J')*A_minus_ind(1)*B_0(:,1);
        PSI_IRS_y(:,ttt+1) = (J*(Q^ttt)*J')*A_minus_ind(2)*B_0(:,2);
        PSI_IRS_f(:,ttt+1) = (J*(Q^ttt)*J')*A_minus_ind(3)*B_0(:,3);
    end
    PSI_0 = [PSI_IRS_m;PSI_IRS_y;PSI_IRS_f];
    PSI_IRS_trim(:,:,kk) = PSI_0;
    
    % Variance Decomposition
    IRF_lim_VD = 100;
    
    PSI_var = zeros(3,3,IRF_lim_VD+1);
    Var_dem_m= zeros(3,IRF_lim_VD+1);
    Var_dem_y= zeros(3,IRF_lim_VD+1);
    Var_dem_f= zeros(3,IRF_lim_VD+1);
    cum_PSI_var = zeros(3,3,IRF_lim_VD+1);
    
    
    for ttt = 0:IRF_lim_VD
        PSI_var(:,:,ttt+1) = ((J*(Q^ttt)*J')* B_0).^2;
        cum_PSI_var(:,:,ttt+1) = sum(PSI_var,3);
        Var_dem_m(:,ttt+1) = (cum_PSI_var(1,:,ttt+1))./sum(cum_PSI_var(1,:,ttt+1));
        Var_dem_y(:,ttt+1) = (cum_PSI_var(2,:,ttt+1))./sum(cum_PSI_var(2,:,ttt+1));
        Var_dem_f(:,ttt+1) = (cum_PSI_var(3,:,ttt+1))./sum(cum_PSI_var(3,:,ttt+1));
    end
    
    
    
    Var_1_dem(:,kk) =[(Var_dem_m(:,1+1)'),(Var_dem_y(:,1+1)'),(Var_dem_f(:,1+1)')]';
    Var_12_dem(:,kk) =[(Var_dem_m(:,12+1)'),(Var_dem_y(:,12+1)'),(Var_dem_f(:,12+1)')]';
    Var_inf_dem(:,kk) =[(Var_dem_m(:,100+1)'),(Var_dem_y(:,100+1)'),(Var_dem_f(:,100+1)')]';
    max_Var_dem(:,kk) =[max(Var_dem_m'),max(Var_dem_y'),max(Var_dem_f')]';
    
    
end

for tttt = 1:IRF_lim+1
    PSI_IRS_m_boot_hi(1,tttt)= prctile(reshape(PSI_IRS_trim(1,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_m_boot_hi(2,tttt)= prctile(reshape(PSI_IRS_trim(2,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_m_boot_hi(3,tttt)= prctile(reshape(PSI_IRS_trim(3,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_y_boot_hi(1,tttt)= prctile(reshape(PSI_IRS_trim(4,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_y_boot_hi(2,tttt)= prctile(reshape(PSI_IRS_trim(5,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_y_boot_hi(3,tttt)= prctile(reshape(PSI_IRS_trim(6,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_f_boot_hi(1,tttt)= prctile(reshape(PSI_IRS_trim(7,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_f_boot_hi(2,tttt)= prctile(reshape(PSI_IRS_trim(8,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    PSI_IRS_f_boot_hi(3,tttt)= prctile(reshape(PSI_IRS_trim(9,tttt,1:ctr),ctr,1),100 - (100-Frac_sol)/2) ;
    
    PSI_IRS_m_boot_lo(1,tttt)= prctile(reshape(PSI_IRS_trim(1,tttt,1:ctr),ctr,1), (100-Frac_sol)/2) ;
    PSI_IRS_m_boot_lo(2,tttt)= prctile(reshape(PSI_IRS_trim(2,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_m_boot_lo(3,tttt)= prctile(reshape(PSI_IRS_trim(3,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_y_boot_lo(1,tttt)= prctile(reshape(PSI_IRS_trim(4,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_y_boot_lo(2,tttt)= prctile(reshape(PSI_IRS_trim(5,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_y_boot_lo(3,tttt)= prctile(reshape(PSI_IRS_trim(6,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_f_boot_lo(1,tttt)= prctile(reshape(PSI_IRS_trim(7,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_f_boot_lo(2,tttt)= prctile(reshape(PSI_IRS_trim(8,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    PSI_IRS_f_boot_lo(3,tttt)= prctile(reshape(PSI_IRS_trim(9,tttt,1:ctr),ctr,1),(100-Frac_sol)/2) ;
    
end

% ================ Variance Decomposition ============================
for kkkk = 1:9
    
    Var_1_decomp_bt(1,kkkk) = prctile(Var_1_dem(kkkk,:),(100-Frac_sol)/2); %lo
    Var_1_decomp_bt(2,kkkk) = prctile(Var_1_dem(kkkk,:),100 - (100-Frac_sol)/2); %hi
    
    Var_12_decomp_bt(1,kkkk) = prctile(Var_12_dem(kkkk,:),(100-Frac_sol)/2); %lo
    Var_12_decomp_bt(2,kkkk) = prctile(Var_12_dem(kkkk,:),100 - (100-Frac_sol)/2); %hi
    
    Var_inf_decomp_bt(1,kkkk) = prctile(Var_inf_dem(kkkk,:),(100-Frac_sol)/2); %lo
    Var_inf_decomp_bt(2,kkkk) = prctile(Var_inf_dem(kkkk,:),100 - (100-Frac_sol)/2); %hi
    
    max_Var_decomp_bt(1,kkkk) = prctile(max_Var_dem(kkkk,:),(100-Frac_sol)/2); %lo
    max_Var_decomp_bt(2,kkkk) = prctile(max_Var_dem(kkkk,:),100 - (100-Frac_sol)/2); %hi
    
end

out_VD.Var_1_decomp_bt = Var_1_decomp_bt;
out_VD.Var_12_decomp_bt = Var_12_decomp_bt;
out_VD.Var_inf_decomp_bt = Var_inf_decomp_bt;
out_VD.max_Var_decomp_bt = max_Var_decomp_bt;

lo = [PSI_IRS_m_boot_lo;PSI_IRS_y_boot_lo;PSI_IRS_f_boot_lo];
hi = [PSI_IRS_m_boot_hi;PSI_IRS_y_boot_hi;PSI_IRS_f_boot_hi];
end