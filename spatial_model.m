%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    Spatial model for multiple regions                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load aggregated Large ensemble
load('Data/TS_mean_reg.mat')
TS_mean_reg = permute(TS_mean_reg,[1 3 2]);

% Load the pre-processed forcings
% Xs: all the covariates in the spatial model
% covariates : 1 vector; f(t); K*2 harmonics; K*2 f(t)*harmonics
load('Data/X_co2.mat')

T=size(TS_mean_reg,2); % 2172 months from 1920 to 2100
R=size(TS_mean_reg,3); % 35 ensemble runs
nreg=size(TS_mean_reg,1); % 47 regions

% Trainning set size chosen: number of realizations = 5
r = 5;

% Creat matrix to store monthly sd, mean trend and residuals 
sd_ind = NaN*ones(nreg,1);
TS_fit_vec = NaN*ones(nreg,T);
res_month_store = NaN * zeros(47,181*r,12);

% ----------Step 1: Estimate the mean trend with linear regression-------%
for reg=1:nreg
    
    
    data=squeeze(TS_mean_reg(reg,:,1:r));
    X=repmat(Xs,r,1);
    % Get parameter estimates
    betahat=(X'*X)\(X'*data(:));
    TS_fit=X*betahat;
    % Store all the residuals
    res_all = data(:) - TS_fit;
    % Estimate the mean trend
    TS_fit_vec(reg,:) = Xs*betahat;
    
    % Store residuals of each month for each region
    for mn = 1:12
        res_80 = [];
        for i = 1:r
            span=mn:12:(12*181);
            index = permute(T*(i-1)+span,[2 1]);
            res = res_all(T*(i-1)+span);
            res_80 = [res_80;res];
        end
    res_month_store(reg,:,mn) = res_80;
    end
end

% ------step 2: Standardize residuals by monthly standard deviation-------%
% Calculate the monthly standard deviation for each region
% and standardize the residuals
std_month = NaN * ones(12,47);
for reg = 1:47
    for mn = 1:12
        data = res_month_store(reg,:,mn);
        std_month(mn,reg) = std(data);
        res_month_store(reg,:,mn) = data/std(data);
    end
end

% Estimate the empirical correlation matrix for standardized residuals
res_stan_LE = reshape(res_month_store,47,[]);
res_stan_LE = permute(res_stan_LE,[2,1]);
cor_month = corr(res_stan_LE);
% save standized residuals as csv file for graphical lasso in R
writematrix(permute(res_stan_LE,[2,1]),'res_stan_LE.csv') 

% -----Step 3: Estimate sparse correlation matrix with Graphical Lasso----%
%---------------------- Uncomment and run in R ------------------------%

% # Download the standardized residuals and perform graphical lass on
% # empirical correlation matrix
% rawdat <- read.csv('res_stan_LE.csv',header=FALSE)
% 
% # Install packages 'spcov' by Bien and Tibshirani (2011)
% install.packages('spcov')
% library(spcov)
% 
% n <- nrow(rawdat)
% p <- ncol(rawdat)
% 
% S <- var(rawdat)
% D <- solve(diag(diag(S)))^(1/2)
% # Sample Correlation Matrix
% R = D
% # compute sparse, positive correlation estimator:
% 
% # Specify the penalty parameter and enable to control the sparsity
% lam <-c(seq(0,1,0.01))
% 
% n_lam <- length(lam)
% # Store the sparse correlation matrix with different penalty
% cov_store <- matrix(NA,47*n_lam,47)
% num_zeros <- matrix(NA,n_lam,1)
% 
% step.size <- 100
% tol <- 1e-3
% P <- matrix(1, p, p)
% diag(P) <- 0
% 
% for (i in 1:n_lam){
%   cat('i=',i,'\n')
%   mm <- spcov(Sigma=R, S=R, lambda=lam[i],
%             step.size=step.size, n.inner.steps=200,
%             thr.inner=0, tol.outer=tol, trace=1)
%   cat('lam=',lam[i],'\n')
%   num_zeros[i,1] = sum(mm$Sigma==0)/(47*47)
%   cov_store[(1+47*(i-1)):(47*i),1:47] = mm$Sigma
% }
% 
% # Write the sparse correlation matrix into csv file
% write.csv(cov_store,'cor_merra2_store_stan.csv')
