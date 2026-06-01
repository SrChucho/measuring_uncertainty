function [B_out] = gen_B_minimal(vec_A, QQ, X_var, X_exp,const_var,data,dates0,lag,ind_spec)

% This is CI Varing B
Vs = data(:,1);
Gold = data(:,2);
T = size(Vs,1);

% Parameter
NQ = size(QQ,3);
eta = (X_var' - repmat(const_var,size(X_var,1),1)'-[vec_A(1:3*lag)';vec_A(3*lag+1:6*lag)';vec_A(6*lag + 1:9*lag)']*X_exp')';
OMEGA= cov(eta);
BB0=chol(OMEGA,'lower'); %Generate lower cholesy factorization of covariance matrix eta

dates3 = 2011 + 7/12;  %2011:07;
dates4 = 2011 + 8/12; %2011:08
dates20 = 1979 + 10/12; %1979:10


%ceiling crisis
index_crsis = [find(dates0 == dates3),find(abs(dates0 -dates4)<10^(-6))];
% 1979
index_79= [find(abs(dates0 -dates20)<10^(-6)):find(abs(dates0 -dates20)<10^(-6))];


% Solve B
parfor iii = 1:NQ
    
    % rotate and normalize B
    B_initial=BB0*QQ(:,:,iii);
    B_est=B_initial*diag(sign(diag(B_initial))); %normalize so that self-irfs are positive
    ehat = eta(:,:)*inv(B_est)';
    
    
    % Compute Correlations
    phi1 = corr(Vs, ehat(:,1));
    phi2 = corr(Vs, ehat(:,3));
    phi1_gold = corr(Gold, ehat(:,1));
    phi2_gold = corr(Gold, ehat(:,3));
    
    
    %============ Event constraint ==============
    
    % gbar4, gbar5 and gbar 6
    ind_ef_event_sign(iii) =  (ehat(index_79,1) >= 0 &ehat(index_79,3) >= 0)& (min(ehat(index_crsis,1)) >= 0&min(ehat(index_crsis,3)) >= 0);
    
    % =========== External Variable Constraint ===========
    ind_ev(iii) =  phi1 <= 0 & phi2 <= 0 & phi1_gold >= 0 & phi2_gold >= 0;
    
    if ind_spec == 1     %cov only
        ind_ef(iii) = 1;
    elseif ind_spec == 2 %cov + external variable
        ind_ef(iii) = ind_ev(iii);
    elseif ind_spec == 3 % cov + external variable + gbar4 - gbar6
        ind_ef(iii) = ind_ev(iii).* ind_ef_event_sign(iii);
    end
    B_store(:,:,iii) = B_est;
    
end

disp('Number of solutions')
disp(sum(ind_ef))

B_out = B_store(:,:,logical(ind_ef));

end