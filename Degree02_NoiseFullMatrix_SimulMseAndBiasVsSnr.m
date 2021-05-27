clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseFullMatrix_SimulMseAndBiasVsSnr'; %Name for the results folder: it should be named after the kind of test performed

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

%% Generate Polynomial and Covariance matrix
% Generate a random covariance
[Sigma,C_atilda,A] = generate_covariance(N,1,'full');
% Generate random roots
r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
% Compute corresponding noise-free polynomial cefficients
a=conj(poly(r)');

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
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
for ii=1:SNR_nsteps
    sigma_a=(sqrt(1/SNRlin(ii)));
    % Compute the expected MSE matrix and bias from the analytic expression
    MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
    MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
    Bias_analytic(:,ii)=bias_analytic(r,a,sigma_a^2*Sigma); % bias (complex augmented)
    Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,ii); % bias (real composite)
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
    Bias_simulated(:,ii)=1/K*sum([err_n(:,:,ii); conj(err_n(:,:,ii))],2);
    Bias_simulated_tilda(:,ii)=1/2*J'*Bias_simulated(:,ii);
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);
D=SNR_nsteps;
for d=1:SNR_nsteps
    subplot(4,D/2,d);
    viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
    plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
    for ii=2:N+1
        plot(real(a_n(ii,:,d)),imag(a_n(ii,:,d)),'.','MarkerSize',1); hold on; % Simulated coeff
    end
%     for ii=1:N
%         ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
%     end
    plot(real(a),imag(a),'*k','MarkerSize',20);
    axis equal;axis(5*[-1,1,-1,1]);
    title("Coefficients");grid on;hold off

    subplot(4,D/2,D+d);
    viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
    plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
    for ii=1:N
        plot(real(r_n(ii,:,d)),imag(r_n(ii,:,d)),'.','MarkerSize',1); hold on; % Simulated roots
    end
    for ii=1:N
    %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
        ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii],d)),[real(r(ii)),imag(r(ii))])
    end
    plot(real(r_mean(:,:,d)),imag(r_mean(:,:,d)),'.b','MarkerSize',15); % Mean of estimated roots
    plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
    axis equal;axis(3*[-1,1,-1,1]);
    title("Roots");grid on;hold off
end

figs(2)=figure(2);

for ii=1:(2*N)
    for jj=ii:(2*N)
        subplot(2*N,2*N,(ii-1)*2*N+jj)
        semilogy(SNR,abs(reshape(MSE_analytic_tilda(ii,jj,:),SNR_nsteps,1)),'-');hold on;
        semilogy(SNR,abs(reshape(MSE_simulated_tilda(ii,jj,:),SNR_nsteps,1)),'o--');
        legend([ strcat("analitic"); strcat("simulated")],'Location','southwest');
        title(strcat("MSE_{",int2str(ii),int2str(jj),"} vs SNR"));grid on;
    end
end

figs(3)=figure(3);

subplot(1,1,1)
leg1=[]; % leg2=[];
for ii=1:(2*N)
    for jj=ii:(2*N)
        semilogy(SNR,abs(reshape(MSE_simulated_tilda(ii,jj,:),SNR_nsteps,1)-reshape(MSE_analytic_tilda(ii,jj,:),SNR_nsteps,1)),'o--'); hold on;
        leg1 =[ leg1; strcat("Error(MSE_{",int2str(ii),int2str(jj),"})")];
    end
end
legend(leg1);
title("MSE error vs SNR");grid on;

figs(4)=figure(4);

for ii=1:N
    subplot(N,2,1+2*(ii-1))
    semilogy(SNR,abs(Bias_analytic_tilda(ii,:)),'-');hold on;
    semilogy(SNR,abs(Bias_simulated_tilda(ii,:)),'o--');
    legend([ strcat("analitic"); strcat("simulated")],'Location','southwest');
    title(strcat("Abs(Bias_{",int2str(ii),"}) vs SNR (real part)"));grid on;
    
    subplot(N,2,2+2*(ii-1))
    semilogy(SNR,abs(Bias_analytic_tilda(N+ii,:)),'-');hold on;
    semilogy(SNR,abs(Bias_simulated_tilda(N+ii,:)),'o--');
    legend([ strcat("analitic"); strcat("simulated")],'Location','southwest');
    title(strcat("Abs(Bias_{",int2str(ii),"}) vs SNR (imaginary part)"));grid on;
end

figs(5)=figure(5);

subplot(1,1,1)
leg1=[]; % leg2=[];
for ii=1:(2*N)
    semilogy(SNR,abs(Bias_analytic_tilda(ii,:)-Bias_simulated_tilda(ii,:)),'o--'); hold on;
    leg1 =[ leg1; strcat("Error(Bias_{",int2str(ii),"})")];
end
legend(leg1);
title("Bias error vs SNR");grid on;

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));