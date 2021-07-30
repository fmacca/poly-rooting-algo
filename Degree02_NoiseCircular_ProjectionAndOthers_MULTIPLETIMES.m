clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_ProjectionAndOthers_MULTIPLETIMES'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
%sigma_a=.01; % variance of the noise on coefficients
SNR = [-12:3:40];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=10;%20; % Number of times we generate a different polynomial

%%
for counter=1:NRUNS
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
    Projection_on_11=zeros(SNR_nsteps,1);
    Norms_r=zeros(SNR_nsteps,1);
    eig_dom=zeros(SNR_nsteps,1);
    Delta_exact=poly2D_discriminant(a);

    % Normality tests
    K_normaltest=2*10^3; % Number of iterations to be used for normality test
%     Gauss_test_Roy=zeros(SNR_nsteps,1); % Matrices to collect the result of Roystest_mod
    Gauss_test_HZ=zeros(SNR_nsteps,1); % Matrices to collect the result of HZmvntest_mod

    % Tests for the mean
%     Ttest_p=zeros(SNR_nsteps,4); % 4 univariate t-tests
    HotT2_p=zeros(SNR_nsteps,1); % Hotelling T^2 test

    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
    %     Bias_analytic(:,ii)=bias_analytic(r(:,d),a(d,:),Sigma); % bias (complex augmented)
    %     Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,d); % bias (real composite)
        % I do things for the projection on orthogonal of [1;1]
    %     Malanobis=inv(MSE_analytic(1:N,1:N,ii)); % This is theMalanobis
    %     metric matrix, for computational reason we do not compute it
        orth=null((MSE_analytic(1:N,1:N,ii)\[1;1])'); % Vector orthodonal to [1;1] in Malanobis metric
        orth_norm=sqrt(orth'*(MSE_analytic(1:N,1:N,ii)\orth)); % Norm of the orthogonal vector
        norm_11=sqrt([1;1]'*(MSE_analytic(1:N,1:N,ii)\[1;1]));
        r_norm=sqrt(r'*(MSE_analytic(1:N,1:N,ii)\r));
        Norms_r(ii)=r_norm;
        Projection(ii)=1/orth_norm*orth'*(MSE_analytic(1:N,1:N,ii)\r);
        Projection_on_11(ii)=1/norm_11*[1;1]'*(MSE_analytic(1:N,1:N,ii)\r);

        eig_dom(ii)=max(abs(sigma_a^2*eig(A'*A)*2));
        for k=1:K
            noise_tilda=sigma_a*A*randn(2*N,1); %Generate colored noise
            a_n(:,k,ii)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
            r_curr=roots(a_n(:,k,ii)); %Compute the roots
            r_n(:,k,ii)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
            err_n(:,k,ii)=r_n(:,k,ii)-r;

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,ii)=1/K*[err_n(:,:,ii); conj(err_n(:,:,ii))]*[err_n(:,:,ii); conj(err_n(:,:,ii))]';
        MSE_simulated_tilda(:,:,ii)=1/4*J'*MSE_simulated(:,:,ii)*J;

%         Gauss_test_Roy(ii)=Roystest_mod([real(r_n(1,1:K_normaltest,ii))' imag(r_n(1,1:K_normaltest,ii))' real(r_n(2,1:K_normaltest,ii))' imag(r_n(2,1:K_normaltest,ii))']);
        Gauss_test_HZ(ii)=HZmvntest_mod([real(r_n(1,:,ii))' imag(r_n(1,:,ii))' real(r_n(2,:,ii))' imag(r_n(2,:,ii))']);
%         [~,Ttest_p(ii,:)]=ttest([real(r_n(1,:,ii))'-real(r(1)) imag(r_n(1,:,ii))'-imag(r(1)) real(r_n(2,:,ii))'-real(r(2)) imag(r_n(2,:,ii))'-imag(r(2))]);
        HotT2_p(ii)=T2Hot1_mod([real(r_n(1,:,ii))'-real(r(1)) imag(r_n(1,:,ii))'-imag(r(1)) real(r_n(2,:,ii))'-real(r(2)) imag(r_n(2,:,ii))'-imag(r(2))]);
    end
    r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

    discr_eig_ratio=abs(Delta_exact)./eig_dom;

    ratio=abs(r(1)-r(2))./abs(Projection);

    close(h); %Close waitbar

%     [SNR' abs(Projection)  abs(r(1)-r(2))./abs(Projection) Gaussianity_test_n Ttest_p];

    % Plots
    figs(1)=figure(1);
    subplot(1,3,1);
    semilogx(abs(Projection),HotT2_p,'x');hold on;grid on;
%     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
    semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'-');
    semilogx([3 3],[0 1],'-');
    xlabel("Abs(projection)");ylabel("pvalue of Hotelling's T^2 test");

    figs(2)=figure(2);
    subplot(1,3,1);
    semilogx(abs(Projection),Gauss_test_HZ,'x');hold on;grid on;
%     semilogx(abs(Projection),Gauss_test_Roy,'x-');
    semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'-');
    semilogx([20 20],[0 1],'-');
    xlabel("Abs(projection)");ylabel("pvalue of mtv. gaussianity test");

    figs(1)=figure(1);
    subplot(1,3,2);
    semilogx(discr_eig_ratio,HotT2_p,'x');hold on;grid on;
%     semilogx(discr_eig_ratio,min(Ttest_p,[],2),'x-');
    semilogx([min(discr_eig_ratio) max(discr_eig_ratio)],[0.05 0.05],'-');
    semilogx([3 3],[0 1],'-');
    xlabel("|\Delta|/\lambda");ylabel("pvalue of Hotelling's T^2 test");

    figs(2)=figure(2);
    subplot(1,3,2);
    semilogx(discr_eig_ratio,Gauss_test_HZ,'x');hold on;grid on;
%     semilogx(discr_eig_ratio,Gauss_test_Roy,'x-');
    semilogx([min(discr_eig_ratio) max(discr_eig_ratio)],[0.05 0.05],'-');
    semilogx([30 30],[0 1],'-');
    xlabel("|\Delta|/\lambda");ylabel("pvalue of mtv. gaussianity test");

    figs(1)=figure(1);
    subplot(1,3,3);
    semilogx(ratio,HotT2_p,'x');hold on;grid on;
%     semilogx(ratio,min(Ttest_p,[],2),'x-');
    semilogx([min(ratio) max(ratio)],[0.05 0.05],'-');
    semilogx([1 1],[0 1],'-');
    xlabel("dist/projection ratio");ylabel("pvalue of Hotelling's T^2 test");

    figs(2)=figure(2);
    subplot(1,3,3);
    semilogx(ratio,Gauss_test_HZ,'x');hold on;grid on;
%     semilogx(ratio,Gauss_test_Roy,'x-');
    semilogx([min(ratio) max(ratio)],[0.05 0.05],'-');
    semilogx([0.2 0.2],[0 1],'-');
    xlabel("dist/projection ratio");ylabel("pvalue of mtv. gaussianity test");

    figs(3)=figure(3);
    subplot(1,2,1);
    loglog(abs(Projection),discr_eig_ratio,'x-');hold on;grid on;
    xlabel("Abs(projection)");ylabel("|\Delta|/\lambda");
    subplot(1,2,2);
    loglog(discr_eig_ratio,ratio,'x-');hold on;grid on;
    xlabel("|\Delta|/\lambda");ylabel("dist/projection ratio");

end

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));