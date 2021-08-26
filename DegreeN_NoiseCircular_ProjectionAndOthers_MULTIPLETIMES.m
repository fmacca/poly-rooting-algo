clear all
close all
clc

addpath('Resources') 
%% Parameters
N=7; % order of the polynomial
%sigma_a=.01; % variance of the noise on coefficients
SNR = 12+[-12:3:40];%[-12:3:40];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^4;%10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
K_normaltest=10^4; % Number of iterations to be used for normality test (since it cannot handle 10^5)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=10; % Number of times we generate a different polynomial

plot_for_thesis=1; % 1 if we want the plots for the publication, 0 if we want the plots for investigation

dataset=[];

%% Generate folder for results
folder_name=strcat('Results/Degree0',int2str(N),'_NoiseCircular_ProjectionAndOthers_MULTIPLETIMES');%'Results/DegreeN_NoiseCircular_ProjectionAndOthers_MULTIPLETIMES'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

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
    n_basis=nchoosek(N,2);
    Projection=zeros(SNR_nsteps,n_basis);
    gamma=zeros(SNR_nsteps,1);

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
    %     Bias_analytic(:,ii)=bias_analytic(r(:,d),a(d,:),Sigma); % bias (complex augmented)
    %     Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,d); % bias (real composite)
        
        % I do things for the projection on orthogonal of the threee spaces
        % corresponding to double roots
    %     Malanobis=inv(MSE_analytic(1:N,1:N,ii)); % This is theMalanobis
    %     metric matrix, for computational reason we do not compute it
        basis=double_root_basis(N);
        for bb=1:n_basis
            orth=null((MSE_analytic(1:N,1:N,ii)\basis(:,:,bb))'); % Vector orthodonal to current basis in Malanobis metric
            orth_norm=sqrt(orth'*(MSE_analytic(1:N,1:N,ii)\orth)); % Norm of the orthogonal vector
            Projection(ii,bb)=1/orth_norm*orth'*(MSE_analytic(1:N,1:N,ii)\r);
        end
        
        gamma(ii)=min(abs(Projection(ii,:)));
        
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

        Gauss_test_HZ(ii)=HZmvntest_mod([real(r_n(:,1:K_normaltest,ii))' imag(r_n(:,1:K_normaltest,ii))']);

        HotT2_p(ii)=T2Hot1_mod([real(r_n(:,:,ii))'-real(r)' imag(r_n(:,:,ii))'-imag(r)']);
    end
    r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
    
    % Save everything into a matrix [counter r(1:N) gamma Gauss_test_HZ
    % HotT2_p]
    dataset=[dataset; [counter*ones(SNR_nsteps,1) ones(SNR_nsteps,1)*conj(r') gamma Gauss_test_HZ HotT2_p]];

    close(h); %Close waitbar

    if plot_for_thesis
        % Plots
        figs(1)=figure(1);
        subplot(1,2,1);
        semilogx(abs(gamma),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
%         semilogx([min(abs(gamma)) max(abs(gamma))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("\gamma(z_0)");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        subplot(1,2,2);
        semilogx(abs(gamma),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
%         semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("\gamma(z_0)");ylabel("p-value of Henze-Zirkler's test");
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
        subplot(1,2,1);
        semilogx(abs(gamma),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
%         semilogx([min(abs(gamma)) max(abs(gamma))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        subplot(1,2,2);
        semilogx(abs(gamma),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
%         semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Henze-Zirkler's test");
        title("Test for normality");
        
        figs(2)=figure(2);
        subplot(1,3,1);
        semilogx(abs(Projection1),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
%         semilogx([min(abs(gamma)) max(abs(gamma))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        subplot(1,3,2);
        semilogx(abs(Projection2),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
%         semilogx([min(abs(gamma)) max(abs(gamma))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        subplot(1,3,3);
        semilogx(abs(Projection3),HotT2_p,'x');hold on;grid on;
    %     semilogx(abs(Projection),min(Ttest_p,[],2),'x-');
%         semilogx([min(abs(gamma)) max(abs(gamma))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([3 3],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Hotelling's T^2 test");
        title("Test for the mean");
        
        figs(3)=figure(3);
        subplot(1,3,1);
        semilogx(abs(Projection1),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
%         semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Henze-Zirkler's test");
        title("Test for normality");
        
        subplot(1,3,2);
        semilogx(abs(Projection2),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
%         semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Henze-Zirkler's test");
        title("Test for normality");
        
        subplot(1,3,3);
        semilogx(abs(Projection3),Gauss_test_HZ,'x');hold on;grid on;
    %     semilogx(abs(Projection),Gauss_test_Roy,'x-');
%         semilogx([min(abs(Projection)) max(abs(Projection))],[0.05 0.05],'r-');
        yline(0.05,'r-');
        semilogx([20 20],[0 1],'b--');
        xlabel("|\gamma(z_0)|");ylabel("p-value of Henze-Zirkler's test");
        title("Test for normality");
    end
end

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));
save(strcat(results_folder,'/dataset'),'dataset');