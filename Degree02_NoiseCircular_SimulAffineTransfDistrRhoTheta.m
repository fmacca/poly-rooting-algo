clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_SimulAffineTransfDistrRhoTheta'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=0.05; % variance of the noise on coefficients
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'circular');
% Generate polynomial
%a=[1.0000 + 0.0000i;   2.5880 - 0.3628i;   1.6431 - 0.4871i];
% a=[1.0000 + 0.0000i;   0;   -1];
% a=[1.0000 + 0.0000i;   0;   -0.5];
%a=[1.0000 + 0.0000i;   0;   -0.001];
% a=[1.0000 + 0.0000i;   0;   -0.25];
a=[1.0000 + 0.0000i;   0;   0];
% Generate random roots
r=roots(a);
[a_transf,~]=poly2D_affinetransf(a);
r_transf=roots(a_transf);

J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K); %Matrix to contain coefficients at every iteration
a_n_transf=zeros(N+1,K); %Matrix to contain transformed coefficients at every iteration
t=zeros(K); %Vector to contain the affine shift at every iteration
% r_n=zeros(N,K); %Matrix to collect roots computed at every iteration
r_n_transf=zeros(N,K); %Matrix to collect roots computed from transformation at every iteration
% r_n_retransf=zeros(N,K); %Matrix to collect roots retrasformed back
% err_n=zeros(N,K); %Matrix to collect the error in roots at every step
for k=1:K
%     noise_tilda=A*randn(2*N,1); %Generate colored noise
    noise_tilda=randn(2*N,1);
%     a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)];
    a_n(:,k)=a+[0;noise_tilda(1)+1i*noise_tilda(N);0];%Add noise to coefficients
    [a_n_transf(:,k),t(k)] = poly2D_affinetransf(a_n(:,k)); %Collect the transformed polynomials
    r_curr=roots(a_n_transf(:,k)); %Compute the roots
    r_n_transf(:,k)=r_curr(order_roots_permutations(r_curr,r_transf)); %Save roots ordered w.r.t. original roots
%     r_n_retransf(:,k)=t(k)+r_n_transf(:,k); %Add the shift to reobtain original roots
%     r_curr=roots(a_n(:,k)); %Compute the roots
%     r_n(:,k)=r_curr(order_roots_permutations(r_curr,r)); %Roots computed without transformation
%     err_n(:,k)=r_n(:,k)-r;

    waitbar(k/K) %Update waitbar
end
% r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
r_transf_mean = mean(r_n_transf,2); %Mean of the roots computed at every iteration

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
axis equal;axis(5*[-1,1,-1,1]);
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
plot(real(a_transf),imag(a_transf),'*k','MarkerSize',20);
axis equal;axis(5*[-1,1,-1,1]);
title("Coefficients of transformed polynomial");grid on;hold off

subplot(2,3,3);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
plot(real(t),imag(t),'.','MarkerSize',1); hold on; % Simulated coeff
% for ii=1:N
%     ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
% end
% plot(real(a),imag(a),'*k','MarkerSize',20);
axis equal;axis(5*[-1,1,-1,1]);
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
plot(real(r_transf_mean),imag(r_transf_mean),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r_transf),imag(r_transf),'*k','MarkerSize',20); % True roots
axis equal;axis(5*[-1,1,-1,1]);
title("Roots of transformed polynomial");grid on;hold off

subplot(2,3,5);
histogram(abs(r_n_transf)); hold on;
plot([1 1],[0 K/10],'r-'); %add red line at rho=1
plot(abs(r_transf),[0 K/10],'g-'); %add green line at rho=rho(r_transf)
title("Distribution of rho");grid on;hold off

subplot(2,3,6);
histogram(angle(r_n_transf(1,:)),'BinEdges',linspace(-pi,pi,100)); hold on;
histogram(angle(r_n_transf(2,:)),'BinEdges',linspace(-pi,pi,100));
plot([angle(r_transf(1)) angle(r_transf(1))],[0 K/40],'g-'); %add green line at theta=theta(r_transf(1))
plot([angle(r_transf(2)) angle(r_transf(2))],[0 K/40],'b-'); %add blue line at theta=theta(r_transf(2))
title("Distribution of theta");grid on;hold off

%% Save workspace and figures to the folder
% savefig(figs,strcat(results_folder,'/figures.fig'))%,'compact');
% clear figs
% save(strcat(results_folder,'/workspace'));