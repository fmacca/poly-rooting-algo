clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_TestManova'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name);
results_folder=strcat(folder_name);

%% Parameters
N=2; % order of the polynomial
%sigma_a=.01; % variance of the noise on coefficients
SNR = [-12:3:40];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
K_normaltest=10^4; % Number of iterations to be used for normality test (since it cannot handle 10^5)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=10;%20; % Number of times we generate a different polynomial

%%
for counter=1:NRUNS
    % Load the dataset with previous simulations (or create a new one)
    if isfile(strcat(results_folder,'/dataset.mat'))
        load(strcat(results_folder,'/dataset.mat'));
    else
         dataset=[];
    end
    
    % Generate Polynomial and Covariance matrix
    % Generate a random circular covariance
    [Sigma,C_atilda,A] = generate_covariance(N,1,'circular');
    % Generate random roots
    r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
    % Compute corresponding noise-free polynomial cefficients
    a=conj(poly(r)');

    % Simulation
    h = waitbar(0,strcat('Simulations in progress ... Please wait... ',int2str(counter),"/",int2str(NRUNS)));
    a_n=zeros(N+1,K,SNR_nsteps); %Matrix to contain coefficients at every iteration
    r_n_sim_analytic=zeros(N,K,SNR_nsteps); %Matrix to collect roots simulated according to analytic formula at every iteration
    r_n=zeros(N,K,SNR_nsteps); %Matrix to collect roots computed at every iteration
    err_n=zeros(N,K,SNR_nsteps); %Matrix to collect the error in roots at every step
    J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
    MSE_analytic=zeros(2*N,2*N,SNR_nsteps);
    MSE_analytic_tilda=zeros(2*N,2*N,SNR_nsteps);
    % Bias_analytic=zeros(2*N,SNR_nsteps);
    % Bias_analytic_tilda=zeros(2*N,SNR_nsteps);
    MSE_simulated=zeros(2*N,2*N,SNR_nsteps);
    MSE_simulated_tilda=zeros(2*N,2*N,SNR_nsteps);
    
    % Proposed indexes of goodness
    Projection=zeros(SNR_nsteps,1);

    % Normality tests
%     Gauss_test_Roy=zeros(SNR_nsteps,1); % Matrices to collect the result of Roystest_mod
    Gauss_test_HZ=3*ones(SNR_nsteps,1); % Matrices to collect the result of HZmvntest_mod
    
    % Test for equality of covariances
    MBox_p=zeros(SNR_nsteps,1);
    
    % Manova test
    Manova_p=zeros(SNR_nsteps,1);
    Manova_d=zeros(SNR_nsteps,1);
    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
        A_MSE=chol(MSE_analytic_tilda(:,:,ii))'; % Matrix to color the syntetic rs
    %     Bias_analytic(:,ii)=bias_analytic(r(:,d),a(d,:),Sigma); % bias (complex augmented)
    %     Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,d); % bias (real composite)
        
        % I do things for the projection on orthogonal of [1;1]
    %     Malanobis=inv(MSE_analytic(1:N,1:N,ii)); % This is theMalanobis
    %     metric matrix, for computational reason we do not compute it
        orth=null((MSE_analytic(1:N,1:N,ii)\[1;1])'); % Vector orthodonal to [1;1] in Malanobis metric
        orth_norm=sqrt(orth'*(MSE_analytic(1:N,1:N,ii)\orth)); % Norm of the orthogonal vector
        
        Projection(ii)=1/orth_norm*orth'*(MSE_analytic(1:N,1:N,ii)\r);

        for k=1:K
            noise_tilda=sigma_a*A*randn(2*N,1); %Generate colored noise
            a_n(:,k,ii)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
            r_curr=roots(a_n(:,k,ii)); %Compute the roots
            r_n(:,k,ii)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
            err_n(:,k,ii)=r_n(:,k,ii)-r;
            
            noise_r_tilda=A_MSE*randn(2*N,1);
            r_n_sim_analytic(:,k,ii)=r+[noise_r_tilda(1:N)+1i*noise_r_tilda(N+1:2*N)];

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,ii)=1/K*[err_n(:,:,ii); conj(err_n(:,:,ii))]*[err_n(:,:,ii); conj(err_n(:,:,ii))]';
        MSE_simulated_tilda(:,:,ii)=1/4*J'*MSE_simulated(:,:,ii)*J;

%         Gauss_test_HZ(ii)=HZmvntest_mod([real(r_n(1,1:K_normaltest,ii))' imag(r_n(1,1:K_normaltest,ii))' real(r_n(2,1:K_normaltest,ii))' imag(r_n(2,1:K_normaltest,ii))']);
        
        MBox_p(ii)=MBoxtest_mod([[ones(K,1) real(r_n(1,:,ii))' imag(r_n(1,:,ii))' real(r_n(2,:,ii))' imag(r_n(2,:,ii))'];[2*ones(K,1) real(r_n_sim_analytic(1,:,ii))' imag(r_n_sim_analytic(1,:,ii))' real(r_n_sim_analytic(2,:,ii))' imag(r_n_sim_analytic(2,:,ii))']]);
        
        [Manova_d(ii),Manova_p(ii)]=manova1([[real(r_n(1,:,ii))' imag(r_n(1,:,ii))' real(r_n(2,:,ii))' imag(r_n(2,:,ii))'];[real(r_n_sim_analytic(1,:,ii))' imag(r_n_sim_analytic(1,:,ii))' real(r_n_sim_analytic(2,:,ii))' imag(r_n_sim_analytic(2,:,ii))']],[ones(K,1);2*ones(K,1)]);
    end
    r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

    % Save everything into a matrix [counter r1 r2 Projection Gauss_test_HZ
    % MBox_p Manova_d Manova_p]
    dataset=[dataset; [counter*ones(SNR_nsteps,1) r(1)*ones(SNR_nsteps,1) r(2)*ones(SNR_nsteps,1) Projection Gauss_test_HZ MBox_p Manova_d Manova_p]];

    % Save the dataset at every iteration
    save(strcat(results_folder,'/dataset'),'dataset');
    
    close(h); %Close waitbar

end

%% Save workspace and figures to the folder
save(strcat(results_folder,'/dataset'),'dataset');
