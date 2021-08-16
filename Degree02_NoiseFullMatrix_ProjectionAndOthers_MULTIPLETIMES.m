clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseFullMatrix_ProjectionAndOthers_MULTIPLETIMES'; %Name for the results folder: it should be named after the kind of test performed

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
K=10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
K_normaltest=10^4; % Number of iterations to be used for normality test (since it cannot handle 10^5)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=10; % Number of times we generate a different polynomial

plot_for_thesis=1; % 1 if we want the plots for the publication, 0 if we want the plots for investigation

dataset=[];

%%
for counter=1:NRUNS
    % Generate Polynomial and Covariance matrix
    % Generate a random circular covariance
    [Sigma,C_atilda,A] = generate_covariance(N,1,'full');
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
    Bias_analytic=zeros(2*N,SNR_nsteps);
    Bias_analytic_tilda=zeros(2*N,SNR_nsteps);
    MSE_simulated=zeros(2*N,2*N,SNR_nsteps);
    MSE_simulated_tilda=zeros(2*N,2*N,SNR_nsteps);
    Bias_simulated=zeros(2*N,SNR_nsteps);
    Bias_simulated_tilda=zeros(2*N,SNR_nsteps);
    
    % Proposed indexes of goodness
    Projection=zeros(SNR_nsteps,1);

    Delta_exact=poly2D_discriminant(a);

    % Normality tests
    Gauss_test_HZ=zeros(SNR_nsteps,1); % Matrices to collect the result of HZmvntest_mod

    % Tests for the mean
    HotT2_p=zeros(SNR_nsteps,1); % Hotelling T^2 test

    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
        Bias_analytic(:,ii)=bias_analytic(r,a,sigma_a^2*Sigma); % bias (complex augmented)
        Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,ii); % bias (real composite)
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

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,ii)=1/K*[err_n(:,:,ii); conj(err_n(:,:,ii))]*[err_n(:,:,ii); conj(err_n(:,:,ii))]';
        MSE_simulated_tilda(:,:,ii)=1/4*J'*MSE_simulated(:,:,ii)*J;
        Bias_simulated(:,ii)=[mean(err_n(:,:,ii),2); conj(mean(err_n(:,:,ii),2))];
        Bias_simulated_tilda(:,ii)=1/2*J'*Bias_simulated(:,ii);
        
        Gauss_test_HZ(ii)=HZmvntest_mod([real(r_n(:,1:K_normaltest,ii))' imag(r_n(:,1:K_normaltest,ii))']);

%         HotT2_p(ii)=T2Hot1_mod([real(r_n(:,:,ii))'-real(r)' imag(r_n(:,:,ii))'-imag(r)']);
        HotT2_p(ii)=T2Hot1_mod([real(err_n(:,:,ii))'-real(Bias_analytic(1:N,ii))' imag(err_n(:,:,ii))'-imag(Bias_analytic(1:N,ii))']);
    end
    r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
    
    % Save everything into a matrix [counter r1 r2 Projection Gauss_test_HZ
    % HotT2_p]
    dataset=[dataset; [counter*ones(SNR_nsteps,1) r(1)*ones(SNR_nsteps,1) r(2)*ones(SNR_nsteps,1) Projection Gauss_test_HZ HotT2_p]];

    close(h); %Close waitbar

%     [SNR' abs(Projection)  abs(r(1)-r(2))./abs(Projection) Gaussianity_test_n Ttest_p];
    if plot_for_thesis
        figs(1)=figure(1);
        subplot(1,2,1);
        semilogx(abs(Projection),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
        semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        subplot(1,2,2);
        semilogx(abs(Projection),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
        semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Henze-Zirkler's test");
        title("Test for normality");
        
        figs(2)=figure(2);
        subplot(3,4,counter);
        viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
        plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
        plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
        axis equal;grid on;

    else % Plots for investigation
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
end

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));
save(strcat(results_folder,'/dataset'),'dataset');