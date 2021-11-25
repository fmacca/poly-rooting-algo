% File adapted to test the idea that things can be done using the
% derivative
clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_SimulDoubleRoot'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=.01; % variance of the noise on coefficients
K=10^4;%10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'circular');
% Generate random double roots
r1=[scale*(2*rand(1,1)-1)+scale*1i*(2*rand(1,1)-1)];
r=[r1; r1];
% Compute corresponding noise-free polynomial cefficients
a=conj(poly(r)');
% Compute polynomial derivative
a_der=1/2*polyder(a);
r_der=roots(a_der);

%% Compute the expected MSE matrix and bias from the analytic expression for the derivative
Sigma_der=Sigma([1 3],[1,3]);
% MSE matrix (complex augmented)
MSE_analytic=mse_analytic(r_der,a_der,(1/4)*Sigma_der);
% MSE matrix (real composite)
J=[eye(N-1) 1i*eye(N-1);eye(N-1) -1i*eye(N-1)]; % Notation change matrix
MSE_analytic_tilda=1/4*J'*MSE_analytic*J;
% % bias (complex augmented)
% Bias_analytic=bias_analytic(r,a,Sigma);
% % bias (real composite)
% Bias_analytic_tilda=1/2*J'*Bias_analytic;

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K); %Matrix to contain coefficients at every iteration
r_n=zeros(N,K); %Matrix to collect roots computed at every iteration
err_n=zeros(N,K); %Matrix to collect the error in roots at every step
Delta_n=zeros(K,1); % Vector to contain the discriminant at every iteration
for k=1:K
    noise_tilda=A*randn(2*N,1); %Generate colored noise
    a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
    r_curr=roots(a_n(:,k)); %Compute the roots
    r_n(:,k)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
    err_n(:,k)=r_n(:,k)-r;
    Delta_n(k)=poly2D_discriminant(a_n(:,k));

    waitbar(k/K) %Update waitbar
end
r_n_unord=[r_n(1,:) r_n(2,:)];
err_n_unord=[err_n(1,:) err_n(2,:)];
J_unord=[eye(N-1) 1i*eye(N-1);eye(N-1) -1i*eye(N-1)]; % Notation change matrix
r_mean = mean(r_n_unord); %Mean of the roots computed at every iteration

MSE_simulated=1/(2*K)*[err_n_unord; conj(err_n_unord)]*[err_n_unord; conj(err_n_unord)]';
MSE_simulated_tilda=1/4*J_unord'*MSE_simulated*J_unord;

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

% subplot(2,2,1);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=2:N+1
    plot(real(a_n(ii,:)),imag(a_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated coeff
end
for ii=1:N
    ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
end
plot(real(a),imag(a),'*k','MarkerSize',20);
axis equal;axis(5*[-1,1,-1,1]);
title("Coefficients");grid on;hold off

figs(2)=figure(2);
% subplot(2,2,2);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n(ii,:)),imag(r_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated roots
end
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
ellipse_plot(0.1*inv(MSE_simulated_tilda),[real(r(1)),imag(r(1))])
ellipse_plot(0.1*inv(MSE_analytic_tilda),[real(r(1)),imag(r(1))])
% end
plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(2*[-1,1,-1,1]);
title("Roots");grid on;hold off

figs(3)=figure(3);
% subplot(2,2,2);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
plot(real(Delta_n),imag(Delta_n),'.','MarkerSize',1); hold on; % Simulated roots
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
% end
% plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
% plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Discriminant");grid on;hold off

figs(4)=figure(4);
% subplot(2,2,2);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
sqrtD=[sqrt(Delta_n);sqrt(Delta_n)*exp(1i*pi)];
plot(real(sqrtD),imag(sqrtD),'.','MarkerSize',1); hold on;
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
% end
% plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
% plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Square root of discriminant");grid on;hold off

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));