clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree03_NoiseCircular_DistanceSNR_OneFixed'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=3; % order of the polynomial
%sigma_a=.01; % variance of the noise on coefficients
SNR = [0:12:36];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
K_normaltest=10^4; % Number of iterations to be used for normality test
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
D=4; % number of shortening steps
% distances=zeros(N,D);
% for ii=1:N
%     distances(ii,:)=abs((4+randn(1))*(1/2).^(0:(D-1)));
% end
distances=abs((4+randn(1))*(1/2).^(0:(D-1)));

confidence=0.05; % Confidence level for normality testing

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,1,'circular');
% Generate a random root baricenter and a random direction for the roots
r1=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
r_fixed=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
dir=rand(1)*2*pi;

%% Simulation
Delta_n=zeros(K,D,SNR_nsteps); % Matrix to contain the discriminant at every iteration
Delta_exact=zeros(D,1); % Matrix to contain the discriminant of exact polynomials
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
r=zeros(N,D); % Exact roots
a=zeros(D,N+1); % Exact coefficients
eig_dom=zeros(1,SNR_nsteps); % Matrix to store dominant eigenvalue of sigma_a^2*Sigma

r_n=zeros(N,K,D,SNR_nsteps); %Matrix to collect roots computed at every iteration

MSE_analytic=zeros(2*N,2*N,D,SNR_nsteps);
MSE_analytic_tilda=zeros(2*N,2*N,D,SNR_nsteps);
% Bias_analytic=zeros(2*N,D,SNR_nsteps);
% Bias_analytic_tilda=zeros(2*N,D,SNR_nsteps);

MSE_simulated=zeros(2*N,2*N,D,SNR_nsteps);
MSE_simulated_tilda=zeros(2*N,2*N,D,SNR_nsteps);

Gaussianity_test_n=zeros(D,SNR_nsteps); % Matrices to collect the result of HZmvntest_mod

for d=1:D
    h = waitbar(0,strcat('Simulations in progress ... Please wait ... ',int2str(d),'/',int2str(D),' ...'));
    r(:,d)=[r_fixed; r1+distances(d)/2*exp(1i*dir); r1-distances(d)/2*exp(1i*dir)];%[r1+distances(:,d)/2*exp(1i*dir); r1-distances(:,d)/2*exp(1i*dir)]; % Set the roots
    a(d,:)=conj(poly(r(:,d))'); % Compute corresponding noise-free polynomial cefficients
    Delta_exact(d)=poly3D_discriminant(a(d,:));
    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute dominant eigenvalue of the covariance matrix
        eig_dom(ii)=max(abs(sigma_a^2*eig(A'*A)*2));
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
            
            Delta_n(k,d,ii)=poly3D_discriminant(a_n);

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,d,ii)=1/K*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]*[err_n(:,:,d,ii); conj(err_n(:,:,d,ii))]';
        MSE_simulated_tilda(:,:,d,ii)=1/4*J'*MSE_simulated(:,:,d,ii)*J;
        
%         Gaussianity_test_n(1,d,ii)=Roystest_mod([real(r_n(1,1:K_normaltest,d,ii))' imag(r_n(1,1:K_normaltest,d,ii))']);
%         Gaussianity_test_n(2,d,ii)=Roystest_mod([real(r_n(2,1:K_normaltest,d,ii))' imag(r_n(2,1:K_normaltest,d,ii))']);
        Gaussianity_test_n(d,ii)=HZmvntest_mod([real(r_n(1,1:K_normaltest,d,ii))' imag(r_n(1,1:K_normaltest,d,ii))' real(r_n(2,1:K_normaltest,d,ii))' imag(r_n(2,1:K_normaltest,d,ii))' real(r_n(3,1:K_normaltest,d,ii))' imag(r_n(3,1:K_normaltest,d,ii))']);
    end
    close(h); %Close waitbar
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
err_mean = mean(err_n,2); %Mean of the error computed at every iteration

discr_eig_ratio=abs(Delta_exact)./eig_dom %An interesting table to look at! discriminant/dominant_eigenvalue

Gaussianity_test_n'

%% Plots
figs(1)=figure(1);
subplot(1,1,1);

viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);

% plot(real(r_mean(:,:,d)),imag(r_mean(:,:,d)),'.b','MarkerSize',15); % Mean of estimated roots

plot([real(r1) real(r(1,:))],[imag(r1) imag(r(1,:))],'-c'); % Direction
plot([real(r1) real(r(2,:))],[imag(r1) imag(r(2,:))],'-c'); % Direction
plot([real(r1) real(r(3,:))],[imag(r1) imag(r(3,:))],'-c'); % Direction
plot(real(r1),imag(r1),'xr','MarkerSize',10); % Baricenter
plot(real(r(1,:)),imag(r(1,:)),'.g','MarkerSize',10); % True roots
plot(real(r(2,:)),imag(r(2,:)),'.r','MarkerSize',10); % True roots
plot(real(r(3,:)),imag(r(3,:)),'.b','MarkerSize',10); % True roots
axis equal;axis([real(r1)-max(max(distances)),real(r1)+max(max(distances)),imag(r1)-max(max(distances)),imag(r1)+max(max(distances))]);
title("Roots");grid on;hold off

%%
figs(2)=figure(2);

for d=1:D
    for ii=1:SNR_nsteps
        subplot(SNR_nsteps,D,d+(ii-1)*D);
        viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
        plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
        
        plot(real(Delta_n(:,d,ii)),imag(Delta_n(:,d,ii)),'.','MarkerSize',1); hold on; % Simulated coeff
        
        color="";% We color with red the title if the P-Value is less than confidence level
        if (Gaussianity_test_n(d,ii)<confidence)
            color="\color{red}";
        end
        axis equal;%axis([min(real(Delta_n(:,d,ii))),max(real(Delta_n(:,d,ii))),min(imag(Delta_n(:,d,ii))),max(imag(Delta_n(:,d,ii)))]);
%         title(strcat(color,"dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii)),"; \Delta/\lambda =",num2str(discr_eig_ratio(d,ii))));grid on;hold off
        title(strcat(color,"dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii))));grid on;hold off
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
        color="";% We color with red the title if the P-Value is less than confidence level
        if (Gaussianity_test_n(d,ii)<confidence)
            color="\color{red}";
        end
        axis equal;
%         title(strcat(color,"dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii)),"; \Delta/\lambda =",num2str(discr_eig_ratio(d,ii))));grid on;hold off
        title(strcat(color,"dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii))));grid on;hold off
    end
end
sgtitle("Roots distributions");

%%
figs(4)=figure(4);

for d=1:D
    for ii=1:SNR_nsteps
        subplot(SNR_nsteps,D,d+(ii-1)*D);
        viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
        plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
        for jj=1:N
            plot(real(err_n(jj,:,d,ii)),imag(err_n(jj,:,d,ii)),'.','MarkerSize',1); hold on; % Simulated roots
        end
        plot(real(err_mean(:,:,d,ii)),imag(err_mean(:,:,d,ii)),'.b','MarkerSize',15); % Mean of estimated roots
        
%         axis equal;axis([min(min(real(r_n(1,:,d,ii))),min(real(r_n(2,:,d,ii)))),max(max(real(r_n(1,:,d,ii))),max(real(r_n(2,:,d,ii)))),min(min(imag(r_n(1,:,d,ii))),min(imag(r_n(2,:,d,ii)))),max(max(imag(r_n(1,:,d,ii))),max(imag(r_n(2,:,d,ii))))]);
        color="";% We color with red the title if the P-Value is less than confidence level
        if (Gaussianity_test_n(d,ii)<confidence)
            color="\color{red}";
        end
        axis equal;
        title(strcat(color,"dist = ",num2str(distances(d)),"; SNR = ",int2str(SNR(ii)),"; \Delta/\lambda =",num2str(discr_eig_ratio(d,ii))));grid on;hold off
    end
end
sgtitle("Error distributions");

%%
figs(5)=figure(5);

leg=[];
for d=1:D
    subplot(1,3,1);
    plot(SNR,abs(reshape(err_mean(1,:,d,:),SNR_nsteps,1)),'o--');hold on;
    leg=[leg; strcat("dist = ",num2str(distances(d)))];
    legend(leg,'Location','northeast');
    title("Average error on root 1 vs SNR");grid on;
    
    subplot(1,3,2);
    plot(SNR,abs(reshape(err_mean(2,:,d,:),SNR_nsteps,1)),'o--');hold on;
    legend(leg,'Location','northeast');
    title("Average error on root 2 vs SNR");grid on;
    
    subplot(1,3,3);
    plot(SNR,abs(reshape(err_mean(3,:,d,:),SNR_nsteps,1)),'o--');hold on;
    legend(leg,'Location','northeast');
    title("Average error on root 3 vs SNR");grid on;
end
for d=1:D
    for ii=1:SNR_nsteps
        subplot(1,3,1);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(1,:,d,ii)),'rx');hold on;
        end
        
        subplot(1,3,2);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(2,:,d,ii)),'rx');hold on;
        end
        
        subplot(1,3,3);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(3,:,d,ii)),'rx');hold on;
        end
    end
end

%%
figs(6)=figure(6);

leg=[];
for d=1:D
    subplot(1,3,1);
    semilogy(SNR,abs(reshape(err_mean(1,:,d,:),SNR_nsteps,1)),'o--');hold on;
    leg=[leg; strcat("dist = ",num2str(distances(d)))];
    legend(leg,'Location','northeast');
    title("Average error on root 1 vs SNR");grid on;
    
    subplot(1,3,2);
    semilogy(SNR,abs(reshape(err_mean(2,:,d,:),SNR_nsteps,1)),'o--');hold on;
    legend(leg,'Location','northeast');
    title("Average error on root 2 vs SNR");grid on;
    
    subplot(1,3,3);
    semilogy(SNR,abs(reshape(err_mean(3,:,d,:),SNR_nsteps,1)),'o--');hold on;
    legend(leg,'Location','northeast');
    title("Average error on root 2 vs SNR");grid on;
end
for d=1:D
    for ii=1:SNR_nsteps
        subplot(1,3,1);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(1,:,d,ii)),'rx');hold on;
        end

        subplot(1,3,2);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(2,:,d,ii)),'rx');hold on;
        end
        
        subplot(1,3,3);
        if (Gaussianity_test_n(d,ii)<confidence)
            plot(SNR(ii),abs(err_mean(3,:,d,ii)),'rx');hold on;
        end
    end
end

%%
figs(7)=figure(7);

for d=1:D
    leg1=[]; leg2=[]; leg3=[];
    for ii=1:SNR_nsteps
        subplot(3,D,d)
        loglog(1:K,abs(cumsum(err_n(1,:,d,ii),2))./repmat(1:K,1,1));
        title(strcat("Root 1, dist = ",num2str(distances(d))));grid on; hold on;
        if d==1
            leg1=[leg1; strcat("SNR = ",int2str(SNR(ii)))];
            legend(leg1,'Location','northeast');
        end
        
        subplot(3,D,D+d)
        loglog(1:K,abs(cumsum(err_n(2,:,d,ii),2))./repmat(1:K,1,1));
        title(strcat("Root 2, dist = ",num2str(distances(d))));grid on; hold on;
        if d==1
            leg2=[leg2; strcat("SNR = ",int2str(SNR(ii)))];
            legend(leg2,'Location','northeast');
        end
        
        subplot(3,D,2*D+d)
        loglog(1:K,abs(cumsum(err_n(3,:,d,ii),2))./repmat(1:K,1,1));
        title(strcat("Root 3, dist = ",num2str(distances(d))));grid on; hold on;
        if d==1
            leg3=[leg3; strcat("SNR = ",int2str(SNR(ii)))];
            legend(leg3,'Location','northeast');
        end
    end
end
sgtitle("Average error at each iteration");

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));