%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    Bias corrected simulations with MERRA2                               %
% Download data from https://disc.gsfc.nasa.gov/datasets?project=MERRA-2  %
% 488 data files from 198001-202008                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Step 1: Apply spatial model-------------------------%

% Trainning set size chosen: number of realizations = 5
r = 5;
% Load aggregated temp data
load('../data/TS_mean_reg.mat')
TS_mean_reg = permute(TS_mean_reg,[1 3 2]);

% Load the forcings(covariates in the spatial model)
load('../data/X_co2.mat')

% Load aggregated reanalysis data MERRA2
load('../data/t2m_merra2.mat')

% Specify the historical period and forecast period for future simulations
time_merra2 = (60*12+9):(60*12+488); % Historical Period 198009 to 202008
time_merra2_pred = (60*12+489):(60*12+488+480); % Forecast period 202009 to 206008
n_pred = length(time_merra2_pred);

res_merra2 = NaN*ones(length(time_merra2),47);
fit_merra2_H = NaN*ones(length(time_merra2),47);
fit_LE_H = NaN*ones(length(time_merra2),47);
res_LE = NaN*ones(length(time_merra2)*r,47);
res_LE_F = NaN*ones(length(time_merra2_pred)*r,47);

fit_merra2_F = NaN*ones(n_pred,47);
fit_LE_F = NaN*ones(length(time_merra2_pred),47);
    
X_co2_H = Xs(time_merra2,:);
X_co2_all = repmat(X_co2_H,r,1);
X_co2_F = Xs(time_merra2_pred,:);
X_co2_F_all = repmat(X_co2_F,r,1);
        
for reg = 1:47
    
    TS_merra2 = t2m_merra2(reg,9:488);
    TS_merra2 = permute(TS_merra2,[2,1]);
    LE = squeeze(TS_mean_reg(reg,time_merra2,1:r));
    LE = LE(:);
    
    LE_F = squeeze(TS_mean_reg(reg,time_merra2_pred,1:r));
    LE_F = LE_F(:);
        
    beta_merra2 = (X_co2_H'*X_co2_H)\(X_co2_H'*TS_merra2);
    beta_LE = (X_co2_all'*X_co2_all)\(X_co2_all'*LE);
    TS_fit_merra2 = X_co2_H*beta_merra2;
    TS_fit_LE = X_co2_all*beta_LE;
    
    TS_err_merra2 = TS_merra2 - TS_fit_merra2;
    TS_err_LE = LE - TS_fit_LE;
        
    res_merra2(:,reg) = TS_err_merra2;
    res_LE(:,reg) = TS_err_LE;
    fit_merra2_H(:,reg) = TS_fit_merra2;
    fit_LE_H(:,reg) = X_co2_H*beta_LE; 
    fit_LE_F(:,reg) = X_co2_F*beta_LE;
    fit_merra2_F(:,reg) = X_co2_F*beta_merra2;
end

% Standardize LE residuals
T = 480;
res_stan_LE = NaN * ones(200/5,12,47);
for reg = 1:47
for mn = 1:12
        res_40 = [];
    for i = 1:r
        span=mn:12:(12*40);
        index = permute(T*(i-1)+span,[2 1]);
        res = res_LE(T*(i-1)+span,reg);
        res_40 = [res_40;res];
    end
    res_stan_LE(:,mn,reg) = res_40; % Historical Period 198009 to 202008
end
end

% Standardize Merra2 residuals
res_stan_merra2 = NaN * ones(40,12,47);
for reg = 1:47 
    for mn = 1:12
        span=mn:12:(12*40);
        index = permute(span,[2 1]);
        res = res_merra2(index,reg);
        res_stan_merra2(:,mn,reg) = res;
    end
end

% Monthly variance for LE and Merra2
var_month_LE = NaN * ones(12,47);
var_month_merra2 = NaN * ones(12,47);
for reg = 1:47
    for mn = 1:12
        data_LE = squeeze(res_stan_LE(:,mn,reg));
        data_merra2 = squeeze(res_stan_merra2(:,mn,reg));
        var_month_LE(mn,reg) = std(data_LE);
        var_month_merra2(mn,reg) = std(data_merra2);
        res_stan_LE(:,mn,reg) = data_LE./std(data_LE);
        res_stan_merra2(:,mn,reg) = data_merra2./std(data_merra2);
    end
end

res_stan_LE = permute(res_stan_LE,[3 1 2]);
res_stan_LE_all = reshape(res_stan_LE,47,[]);
cor_month_LE = corr(permute(res_stan_LE_all,[2,1]));

res_stan_merra2 = permute(res_stan_merra2,[3 1 2]);
res_stan_merra2_all = reshape(res_stan_merra2,47,[]);
cor_month_merra2 = corr(permute(res_stan_merra2_all,[2,1]));

% writematrix(permute(res_stan_LE_all,[2,1]),'res_stan_LE.csv') 
% writematrix(permute(res_stan_merra2_all,[2,1]),'res_stan_merra2.csv') 

% -------------Step 2:Generate bias corrected simulations-----------------%

% Load sparse correlation matrix for LE and MERRA2
merra2_cor_store = csvread('../Rdata/cor_merra2_store_stan.csv');
LE_cor_all_store = csvread('../Rdata/cor_LE_store_stan.csv');

lam_merra2 = (0.1:0.1:1); % Penaly used for MERRA2
lam_LE = (0:0.01:1); % Penaly used for LE
N_lam_merra2 = length(lam_merra2);
N_lam_LE = length(lam_LE);
result_merra2_cor = cell(N_lam_merra2, 1);
result_LE_all_cor = cell(N_lam_LE, 1);

for i = 1:N_lam_merra2
    new_cor_merra2 = merra2_cor_store((1+47*(i-1)):(47*i),1:47);
    result_merra2_cor{i} = new_cor_merra2; 
end

for i = 1:N_lam_LE
    new_cor_LE = LE_cor_all_store((1+47*(i-1)):(47*i),1:47);
    result_LE_all_cor{i} = new_cor_LE; 
end

% Calculating how many zero entries in the sparse correlation matrix
epsilon = 1e-20;
zero_new_merra2 = zeros(N_lam_merra2,1);
zero_new_LE_all = zeros(N_lam_LE,1);

for i=1:N_lam_merra2
    W_merra2 = result_merra2_cor{i};
    zero_new_merra2(i) = sum(abs(W_merra2(:))<epsilon)/(47*47);
end

for i=1:N_lam_LE
    W_LE = result_LE_all_cor{i};
    zero_new_LE_all(i) = sum(abs(W_LE(:))<epsilon)/(47*47);
end

% Choose lasso penalty lambda for achieving 60% sparsity
i_merra2 = 2;
lam_merra2(i_merra2
zero_new_merra2(i_merra2);

i_LE = 16;
lam_LE(i_LE);
zero_new_LE_all(i_LE);

% 60% sparse correlation matrix estimated for LE and MERRA2
Sigma_merra2 = result_merra2_cor{i_merra2};
Sigma_LE_all = result_LE_all_cor{i_LE};
mean_MLR = zeros(47,1);
fit_LE_F = fit_LE_F';

% Generate 30 bias corrected simulations
nsim = 30;
merra2_F_store = NaN*ones(nreg,length(time_merra2_pred),nsim);
merra2_F_store_LE = NaN*ones(nreg,length(time_merra2_pred),nsim);

for i = 1:nsim
    
    for mn = 1:12
        cov_LE_month = NaN * ones(47,47);
        cov_merra2_month = NaN * ones(47,47);
        for reg1 = 1:47
            for reg2 = 1:47
                cov_LE_month(reg1,reg2) = Sigma_LE_all(reg1,reg2) * var_month_LE(mn,reg1) * var_month_LE(mn,reg2);
                cov_merra2_month(reg1,reg2) = Sigma_merra2(reg1,reg2) * var_month_merra2(mn,reg1) * var_month_merra2(mn,reg2);
            end
        end
        TS_err_MLR_month = mvnrnd(mean_MLR, cov_LE_month, length(time_merra2_pred)/12);
        span=mn:12:(12*81);
    
        LE_sim_F = fit_LE_F(span,:) + TS_err_MLR_month;
    
        W_LE_F_month = LE_sim_F;
        fit_LE_F_month = fit_LE_F(span,:);
        fit_merra2_F_month = fit_merra2_F(span,:);
    
        % Generate simulations with bias correction formula
        merra2_F = fit_merra2_F_month'+ chol(cov_merra2_month)*(inv(chol(cov_LE_month)))'*(W_LE_F_month' - fit_LE_F_month');
        merra2_F_store(:,span,i) = merra2_F;
     end
    
end
