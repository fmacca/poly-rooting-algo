clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_SimulAffineTransf'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=.01; % variance of the noise on coefficients
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'circular');
% Generate "roots of 1" polynomial
a=[1; 0; 1];
% Generate random roots
r=roots(a);

J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K); %Matrix to contain coefficients at every iteration
a_n_transf=zeros(N+1,K); %Matrix to contain transformed coefficients at every iteration
t=zeros(K); %Vector to contain the affine shift at every iteration
r_n=zeros(N,K); %Matrix to collect roots computed at every iteration
r_n_transf=zeros(N,K); %Matrix to collect roots computed from transformation at every iteration
err_n=zeros(N,K); %Matrix to collect the error in roots at every step
for k=1:K
    noise_tilda=A*randn(2*N,1); %Generate colored noise
    a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
    [a_n_transf(:,k),t(k)] = poly2D_affinetransf(a_n(:,k)); %Collect the transformed polynomials
    r_curr=roots(a_n_transf(:,k)); %Compute the roots
    r_n_transf(:,k)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
    r_n(:,k)=t(k)+r_n_transf(:,k); %Add the shift to reobtain original roots
    err_n(:,k)=r_n(:,k)-r;

    waitbar(k/K) %Update waitbar
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

subplot(2,3,1);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=2:N+1
    plot(real(a_n(ii,:)),imag(a_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated coeff
end
for ii=1:N
    ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
end
plot(real(a),imag(a),'*k','MarkerSize',20);
axis equal;axis(1.5*[-1,1,-1,1]);
title("Coefficients of original polynomial");grid on;hold off

subplot(2,3,2);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=2:N+1
    plot(real(a_n_transf(ii,:)),imag(a_n_transf(ii,:)),'.','MarkerSize',1); hold on; % Simulated coeff
end
% for ii=1:N
%     ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
% end
plot(real(a),imag(a),'*k','MarkerSize',20);
axis equal;axis(1.5*[-1,1,-1,1]);
title("Coefficients of transformed polynomial");grid on;hold off

subplot(2,3,3);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
plot(real(t),imag(t),'.','MarkerSize',1); hold on; % Simulated coeff
% for ii=1:N
%     ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
% end
% plot(real(a),imag(a),'*k','MarkerSize',20);
axis equal;axis(1.5*[-1,1,-1,1]);
title("Shift paramenter t");grid on;hold off

subplot(2,3,4);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n_transf(ii,:)),imag(r_n_transf(ii,:)),'.','MarkerSize',1); hold on; % Simulated roots
end
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
% end
plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Roots of transformed polynomial");grid on;hold off

subplot(2,3,5);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n(ii,:)),imag(r_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated roots
end
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
% end
plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(1.5*[-1,1,-1,1]);
title("Roots of the original polynomial");grid on;hold off

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));