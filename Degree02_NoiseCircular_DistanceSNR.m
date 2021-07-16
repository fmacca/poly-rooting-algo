clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_DistanceSNR'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
%sigma_a=.01; % variance of the noise on coefficients
SNR = [0:12:36];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
D=4; % number of shortening steps
distances=abs((4+randn(1))*(1/2).^(0:(D-1)));

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,1,'circular');
% Generate a random root baricenter and a random direction for the roots
r1=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
dir=rand(1)*2*pi;

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
Delta_n=zeros(K,D,SNR_nsteps); % Matrix to contain the discriminant at every iteration
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
r=zeros(N,D); % Exact roots
a=zeros(D,N+1); % Exact coefficients

r_n=zeros(N,K,D,SNR_nsteps); %Matrix to collect roots computed at every iteration

MSE_analytic=zeros(2*N,2*N,D,SNR_nsteps);
MSE_analytic_tilda=zeros(2*N,2*N,D,SNR_nsteps);
% Bias_analytic=zeros(2*N,D,SNR_nsteps);
% Bias_analytic_tilda=zeros(2*N,D,SNR_nsteps);

MSE_simulated=zeros(2*N,2*N,D,SNR_nsteps);
MSE_simulated_tilda=zeros(2*N,2*N,D,SNR_nsteps);

for d=1:D
    r(:,d)=[r1+distances(d)/2*exp(1i*dir); r1-distances(d)/2*exp(1i*dir)]; % Set the two roots
    a(d,:)=conj(poly(r(:,d))'); % Compute corresponding noise-free polynomial cefficients
    
    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,d,ii)=mse_analytic(r(:,d),a(d,:),sigma_a^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,d,ii)=1/4*J'*MSE_analytic(:,:,d,ii)*J; % MSE matrix (real composite)
%         Bias_analytic(:,d)=bias_analytic(r(:,d),a(d,:),Sigma); % bias (complex augmented)
%         Bias_analytic_tilda(:,d)=1/2*J'*Bias_analytic(:,d); % bias (real composite)
        for k=1:K
            noise_tilda=sigma_a*A*randn(2*N,1); %Generate colored noise
            a_n=conj(a(d,:)')+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
            r_curr=roots(a_n); %Compute the roots
            r_n(:,k,d,ii)=r_curr(order_roots_permutations(r_curr,r(:,d))); %Save roots ordered w.r.t. original roots
            err_n(:,k,d,ii)=r_n(:,k,d,ii)-r(:,d);
            
            Delta_n(k,d,ii)=poly2D_discriminant(a_n);

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,d,ii)=1/K*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]';
        MSE_simulated_tilda(:,:,d,ii)=1/4*J'*MSE_simulated(:,:,d,ii)*J;
    end
    
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
close(h); %Close waitbar

%% Plots
figs(1)=figure(1);
subplot(1,1,1);

viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);

% plot(real(r_mean(:,:,d)),imag(r_mean(:,:,d)),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r(1,:)),imag(r(1,:)),'.b','MarkerSize',10); % True roots
plot(real(r(2,:)),imag(r(2,:)),'.r','MarkerSize',10); % True roots
axis equal;axis([real(r1)-max(distances),real(r1)+max(distances),imag(r1)-max(distances),imag(r1)+max(distances)]);
title("Roots");grid on;hold off

%%
figs(2)=figure(2);

for d=1:D
    for ii=1:SNR_nsteps
        subplot(SNR_nsteps,D,d+(ii-1)*D);
        viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
        plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
        
        plot(real(Delta_n(:,d,ii)),imag(Delta_n(:,d,ii)),'.','MarkerSize',1); hold on; % Simulated coeff
                
        axis equal;%axis([min(real(Delta_n(:,d,ii))),max(real(Delta_n(:,d,ii))),min(imag(Delta_n(:,d,ii))),max(imag(Delta_n(:,d,ii)))]);
        title(strcat("dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii))));grid on;hold off
    end
end
sgtitle("Discriminant distribution");

%%
figs(3)=figure(3);

for d=1:D
    for ii=1:SNR_nsteps
        subplot(SNR_nsteps,D,d+(ii-1)*D);
        viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
        plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
        for jj=1:N
            plot(real(r_n(jj,:,d,ii)),imag(r_n(jj,:,d,ii)),'.','MarkerSize',1); hold on; % Simulated roots
        end
        for jj=1:N
        %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
            ellipse_plot(0.1*inv(MSE_analytic_tilda([jj N+jj],[jj N+jj],d,ii)),[real(r(jj,d)),imag(r(jj,d))])
        end
        plot(real(r_mean(:,:,d,ii)),imag(r_mean(:,:,d,ii)),'.b','MarkerSize',15); % Mean of estimated roots
        plot(real(r(:,d)),imag(r(:,d)),'*k','MarkerSize',20); % True roots
%         axis equal;axis([min(min(real(r_n(1,:,d,ii))),min(real(r_n(2,:,d,ii)))),max(max(real(r_n(1,:,d,ii))),max(real(r_n(2,:,d,ii)))),min(min(imag(r_n(1,:,d,ii))),min(imag(r_n(2,:,d,ii)))),max(max(imag(r_n(1,:,d,ii))),max(imag(r_n(2,:,d,ii))))]);
        axis equal;
        title(strcat("dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii))));grid on;hold off
    end
end
sgtitle("Roots distributions");


%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));