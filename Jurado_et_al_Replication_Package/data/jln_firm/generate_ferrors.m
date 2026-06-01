% -------------------------------------------------------------------------
% Generate forecast errors 
% -------------------------------------------------------------------------

% Load data
clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
xt          = data;

% Estimate macro factors
[e,fhat,lf,vf] = factors(xt,20,2,2);
[e,ghat,lg,vg] = factors(xt.^2,20,2,2);
ft = [fhat,fhat(:,1).^2,ghat(:,1)]; %predictor set

% Firm data
[data,txt] = xlsread('firmpanel.xlsx');
names      = txt;
gvkey      = unique(data(:,1));
dates      = unique(data(:,2));
T          = length(dates);
N          = length(gvkey);
dpretaxy   = reshape(data(:,3),T,N);
dpretax    = reshape(data(:,4),T,N);
gassets    = reshape(data(:,5),T,N);
gsales     = reshape(data(:,6),T,N);
gprice     = reshape(data(:,7),T,N);
xt         = dpretaxy;

% Take quarterly averages of macro predictors
mdates(1) = 1960.25;
for i = 1:3:length(ft)-1
    fm(1+(i-1)/3,:) = mean(ft(i:i+2,:));
    if i>1; mdates(1+(i-1)/3)  = mdates((i-1)/3)+0.25; end;
end
a  = find(mdates==dates(1));
b  = find(mdates==dates(end));
fm = fm(a:b,:);

% Estimate firm factors
[e,fhat,lf,vf] = factors(xt,20,2,2);
[e,ghat,lg,vg] = factors(xt.^2,20,2,2);
outf     = 'mR2_fhat.out';
outg     = 'mR2_ghat.out';
[R2,mR2] = mrsq(fhat,lf,vf,gvkey,outf);
[R2,mR2] = mrsq(ghat,lg,vg,gvkey,outg);
ft       = [fhat,fhat.^2,ghat,fm]; %predictor set

% Generate forecast errors for yt
yt     = zscore(xt);
[T,N]  = size(yt);
py     = 2;
pz     = 1;
p      = max(py,pz);
q      = fix(4*(T/100)^(2/9));
ybetas = zeros(1+py+pz*size(ft,2),N);
for i = 1:N
    X    = [ones(T,1),mlags(yt(:,i),py),mlags(ft,pz)];
    reg  = nwest(yt(p+1:end,i),X(p+1:end,:),q);
    pass = abs(reg.tstat(py+2:end)) > 2.575; % hard threshold
    keep = [ones(1,py+1)==1,pass'];
    Xnew = X(:,keep);
    reg  = nwest(yt(p+1:end,i),Xnew(p+1:end,:),q);
    vyt(:,i)       = reg.resid; % forecast errors
    ybetas(keep,i) = reg.beta;   
    fmodels(:,i)   = pass; %chosen predictors
end

% Generate AR(4) errors for ft
[T,R]  = size(ft);
pf     = 2;
q      = fix(4*(T/100)^(2/9));
fbetas = zeros(R,pf+1);
for i = 1:R
   X   = [ones(T,1),mlags(ft(:,i),pf)];
   reg = nwest(ft(pf+1:end,i),X(pf+1:end,:),q);
   vft(:,i)    = reg.resid;
   fbetas(i,:) = reg.beta';
end

% Save data
[T,N]  = size(vyt);
ybetas = ybetas';
dates  = 1900+(59+1/4:1/4:111.25)';
dates  = dates(end-T+1:end);
save ferrors dates vyt vft names gvkey ybetas fbetas py pz pf ft xt fmodels

% Also write to .txt file for R code
dlmwrite('vyt.txt',vyt,'delimiter','\t','precision',17);
dlmwrite('vft.txt',vft,'delimiter','\t','precision',17);