clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_SimulBiasVsIterDecreasingDistance'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=.002; % variance of the noise on coefficients
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
D=10; % number of shortening steps
distances=(1/2).^(0:(D-1));

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'circular');
% Generate a random root and a random direction for the second root
r1=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
dir=rand(1)*2*pi;
% % Generate random roots
% r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
% % Compute corresponding noise-free polynomial cefficients
% a=conj(poly(r)');

%% Compute the expected MSE matrix and bias from the analytic expression
% % MSE matrix (complex augmented)
% MSE_analytic=mse_analytic(r,a,Sigma);
% % MSE matrix (real composite)
% J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
% MSE_analytic_tilda=1/4*J'*MSE_analytic*J;
% % bias (complex augmented)
% Bias_analytic=bias_analytic(r,a,Sigma);
% % bias (real composite)
% Bias_analytic_tilda=1/2*J'*Bias_analytic;

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K,D); %Matrix to contain coefficients at every iteration
r_n=zeros(N,K,D); %Matrix to collect roots computed at every iteration
err_n=zeros(N,K,D); %Matrix to collect the error in roots at every step
for d=1:D
    r=[r1; r1+distances(d)*exp(1i*dir)]; % Set the two roots
    a=conj(poly(r)'); % Compute corresponding noise-free polynomial cefficients
    for k=1:K
        noise_tilda=A*randn(2*N,1); %Generate colored noise
        a_n(:,k,d)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
        r_curr=roots(a_n(:,k,d)); %Compute the roots
        r_n(:,k,d)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
        err_n(:,k,d)=r_n(:,k,d)-r;

        waitbar((d*K+k)/(K*D)) %Update waitbar
    end
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

for d=1:D
    subplot(4,D/2,d);
    viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
    plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
    for ii=2:N+1
        plot(real(a_n(ii,:,d)),imag(a_n(ii,:,d)),'.','MarkerSize',1); hold on; % Simulated coeff
    end
    for ii=1:N
        ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
    end
%     plot(real(a),imag(a),'*k','MarkerSize',20);
    axis equal;axis(5*[-1,1,-1,1]);
    title("Coefficients");grid on;hold off

    subplot(4,D/2,D+d);
    viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
    plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
    for ii=1:N
        plot(real(r_n(ii,:,d)),imag(r_n(ii,:,d)),'.','MarkerSize',1); hold on; % Simulated roots
    end
%     for ii=1:N
%     %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%         ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
%     end
    plot(real(r_mean(:,:,d)),imag(r_mean(:,:,d)),'.b','MarkerSize',15); % Mean of estimated roots
%     plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
    axis equal;axis(3*[-1,1,-1,1]);
    title("Roots");grid on;hold off
end

figs(2)=figure(2);

subplot(1,1,1)
leg1=[]; % leg2=[];
for d=1:D
    loglog(1:K,abs(cumsum(err_n(:,:,d),2))./[1:K;1:K]);
    title("Average error (cumulative) vs iteration");grid on; hold on;
    
    % Organizing Legend
    for ii=1:N
        leg1 =[ leg1; strcat("Avg err root ",int2str(ii),"; dist = 1/2^",int2str(d-1))];
    %     leg2 = [leg2; strcat("An. bias root ",int2str(ii))];
    end
end
legend(leg1);

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));