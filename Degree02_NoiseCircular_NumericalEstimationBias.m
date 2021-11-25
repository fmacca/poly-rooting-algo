clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_NumericalEstimationBias'; %Name for the results folder: it should be named after the kind of test performed

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
K=10^4;%10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
% K_normaltest=2*10^3; % Number of iterations to be used for normality test
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
D=10; % number of shortening steps
distances=abs(6*(1/2).^(0:(D-1))); %distances=abs((4+randn(1))*(1/2).^(0:(D-1)));

% confidence=0.05; % Confidence level for normality testing

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,1,'circular');
% % Generate a random root baricenter and a random direction for the roots
% r1=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
% dir=rand(1)*2*pi;

% Select a root baricenter and a direction for the roots
% r1=1*exp(1i*pi/4);
% dir=3/4*pi;

% r1=1*exp(1i*pi/4);
% dir=1/4*pi;

r1=0;
dir=1/4*pi;

%% Simulation
% Delta_n=zeros(K,D,SNR_nsteps); % Matrix to contain the discriminant at every iteration
% Delta_exact=zeros(D,1); % Matrix to contain the discriminant of exact polynomials
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
r=zeros(N,D); % Exact roots
a=zeros(D,N+1); % Exact coefficients
% eig_dom=zeros(1,SNR_nsteps); % Matrix to store dominant eigenvalue of sigma_a^2*Sigma

r_n=zeros(N,K,D,SNR_nsteps); %Matrix to collect roots computed at every iteration
err_n=zeros(N,K,D,SNR_nsteps);

MSE_analytic=zeros(2*N,2*N,D,SNR_nsteps);
MSE_analytic_tilda=zeros(2*N,2*N,D,SNR_nsteps);
% Bias_analytic=zeros(2*N,D,SNR_nsteps);
% Bias_analytic_tilda=zeros(2*N,D,SNR_nsteps);

MSE_simulated=zeros(2*N,2*N,D,SNR_nsteps);
MSE_simulated_tilda=zeros(2*N,2*N,D,SNR_nsteps);

% Gaussianity_test_n=zeros(D,SNR_nsteps); % Matrices to collect the result of HZmvntest_mod

for d=1:D
    h = waitbar(0,strcat('Simulations in progress ... Please wait ... ',int2str(d),'/',int2str(D),' ...'));
    r(:,d)=[r1+distances(d)/2*exp(1i*dir); r1-distances(d)/2*exp(1i*dir)]; % Set the two roots
    a(d,:)=conj(poly(r(:,d))'); % Compute corresponding noise-free polynomial cefficients
%     Delta_exact(d)=poly2D_discriminant(a(d,:));
    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute dominant eigenvalue of the covariance matrix
%         eig_dom(ii)=max(abs(sigma_a^2*eig(A'*A)*2));
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
            err_n_phase(:,k,d,ii)=angle(r_n(:,k,d,ii)./r(:,d)); %phase shift
            
%             Delta_n(k,d,ii)=poly2D_discriminant(a_n);

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,d,ii)=1/K*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]';
        MSE_simulated_tilda(:,:,d,ii)=1/4*J'*MSE_simulated(:,:,d,ii)*J;
        
%         Gaussianity_test_n(1,d,ii)=Roystest_mod([real(r_n(1,1:K_normaltest,d,ii))' imag(r_n(1,1:K_normaltest,d,ii))']);
%         Gaussianity_test_n(2,d,ii)=Roystest_mod([real(r_n(2,1:K_normaltest,d,ii))' imag(r_n(2,1:K_normaltest,d,ii))']);
%         Gaussianity_test_n(d,ii)=HZmvntest_mod([real(r_n(1,:,d,ii))' imag(r_n(1,:,d,ii))' real(r_n(2,:,d,ii))' imag(r_n(2,:,d,ii))']);
    end
    close(h); %Close waitbar
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
err_mean = mean(err_n,2); %Mean of the error computed at every iteration
err_phase_mean = mean(err_n_phase,2); %Mean of the phase error computed at every iteration

% discr_eig_ratio=abs(Delta_exact)./eig_dom %An interesting table to look at! discriminant/dominant_eigenvalue

% Gaussianity_test_n'

%% Plot
figs(1)=figure(1);
subplot(1,1,1);

viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);

% plot(real(r_mean(:,:,d)),imag(r_mean(:,:,d)),'.b','MarkerSize',15); % Mean of estimated roots

plot([real(r(1,:)) real(r(2,:))],[imag(r(1,:)) imag(r(2,:))],'-c'); % Direction
plot(real(r1),imag(r1),'xr','MarkerSize',10); % Baricenter
plot(real(r(1,:)),imag(r(1,:)),'.b','MarkerSize',10); % True roots
plot(real(r(2,:)),imag(r(2,:)),'.r','MarkerSize',10); % True roots
plot(real(r(2,:)),imag(r(2,:)),'.r','MarkerSize',10); % True roots
axis equal;axis([real(r1)-max(distances),real(r1)+max(distances),imag(r1)-max(distances),imag(r1)+max(distances)]);
title("Roots");grid on;hold off

%%
figs(2)=figure(2);

for n=1:N
    subplot(1,N,n);
    leg=[];
    for ii=1:SNR_nsteps
        to_plot=err_mean(n,1,:,ii);
        to_plot=permute(to_plot,[1 3 2]);
        plot(distances,abs(to_plot),'--x');
        xlabel("Distance");
        leg=[leg; strcat("SNR = ",num2str(SNR(ii)),"dB")];
        hold on;grid on;
    end
    title(strcat("Root ",num2str(n)));grid on;hold off
end
legend(leg,'Location','northeast');
sgtitle("Bias Absolute Value");

% %%
% figs(3)=figure(3);
% 
% for n=1:N
%     subplot(1,N,n);
%     leg=[];
%     for ii=1:SNR_nsteps
%         to_plot=err_mean(n,1,:,ii);
%         to_plot=permute(to_plot,[1 3 2]);
%         plot(distances,angle(to_plot),'--x');
%         xlabel("Distance");
%         leg=[leg; strcat("SNR = ",num2str(SNR(ii)),"dB")];
%         hold on;grid on;
%     end
%     title(strcat("Root ",num2str(n)));grid on;hold off
% end
% legend(leg,'Location','northeast');
% sgtitle("Bias Phase");

%%
figs(3)=figure(3);

for n=1:N
    subplot(1,N,n);
    leg=[];
    for ii=1:SNR_nsteps
        to_plot=err_phase_mean(n,1,:,ii);
        to_plot=permute(to_plot,[1 3 2]);
        plot(distances,to_plot,'--x');
        xlabel("Distance");
        leg=[leg; strcat("SNR = ",num2str(SNR(ii)),"dB")];
        hold on;grid on;
    end
    title(strcat("Root ",num2str(n)));grid on;hold off
end
legend(leg,'Location','northeast');
sgtitle("Bias Phase");

%%
figs(4)=figure();

for n=1:N
    subplot(1,N,n);
    leg=[];
    for ii=1:SNR_nsteps
        to_plot=err_mean(n,1,:,ii);
        to_plot=permute(to_plot,[1 3 2]);
        plot(distances,real(to_plot*conj(exp(1i*dir))),'--x');
        xlabel("Distance");
        leg=[leg; strcat("SNR = ",num2str(SNR(ii)),"dB")];
        hold on;grid on;
    end
    title(strcat("Root ",num2str(n)));grid on;hold off
end
legend(leg,'Location','northeast');
sgtitle("Projection of bias on the direction indentified by the roots");

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));