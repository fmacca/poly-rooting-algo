clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseFullMatrix_SimulBiasVsIterOne'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);

%% Parameters
N=2; % order of the polynomial
sigma_a=.01; % variance of the noise on coefficients
K=10^5; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square

%% Generate Polynomial and Covariance matrix
% Generate a random circular covariance
[Sigma,C_atilda,A] = generate_covariance(N,sigma_a,'full');
% Generate random roots
r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
% Compute corresponding noise-free polynomial cefficients
a=conj(poly(r)');

%% Compute the expected MSE matrix and bias from the analytic expression
% MSE matrix (complex augmented)
MSE_analytic=mse_analytic(r,a,Sigma);
% MSE matrix (real composite)
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
MSE_analytic_tilda=1/4*J'*MSE_analytic*J;
% bias (complex augmented)
Bias_analytic=bias_analytic(r,a,Sigma);
% bias (real composite)
Bias_analytic_tilda=1/2*J'*Bias_analytic;

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
a_n=zeros(N+1,K); %Matrix to contain coefficients at every iteration
r_n=zeros(N,K); %Matrix to collect roots computed at every iteration
err_n=zeros(N,K); %Matrix to collect the error in roots at every step
for k=1:K
    noise_tilda=A*randn(2*N,1); %Generate colored noise
    a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
    r_curr=roots(a_n(:,k)); %Compute the roots
    r_n(:,k)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
    err_n(:,k)=r_n(:,k)-r;

    waitbar(k/K) %Update waitbar
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

subplot(2,2,1);
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

subplot(2,2,2);
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n(ii,:)),imag(r_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated roots
end
for ii=1:N
%     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
    ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
end
plot(real(r_mean),imag(r_mean),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(3*[-1,1,-1,1]);
title("Roots");grid on;hold off

% subplot(2,2,3)
% loglog(1:K,abs(real(cumsum(err_n,2)))./[1:K;1:K]); hold on
% title("Abs(Re(Average error)) vs iteration");grid on;
% loglog(1:K,repmat(abs(real(Bias_analytic(1:N))),1,K));
% % Organizing Legend
% leg1=[]; leg3=[];
% for ii=1:N
%     leg1 =[ leg1; strcat("Avg err r_",int2str(ii),", real part")];
%     leg3 = [leg3; strcat("An. bias r_",int2str(ii),", real part")];
% end
% legend([leg1; leg3],"Location","Southwest");
% 
% subplot(2,2,4)
% loglog(1:K,abs(imag(cumsum(err_n,2)))./[1:K;1:K]); hold on
% title("Abs(Im(Average error)) vs iteration");grid on;
% loglog(1:K,repmat(abs(imag(Bias_analytic(1:N))),1,K));
% % Organizing Legend
% leg2=[]; leg4=[];
% for ii=1:N
%     leg2 = [leg2; strcat("Avg err r_",int2str(ii),", imag part")];
%     leg4 = [leg4; strcat("An. bias r_",int2str(ii),", imag part")];
% end
% legend([leg2; leg4],"Location","Southwest");

subplot(2,2,3)
loglog(1:K,abs(cumsum(err_n,2))./[1:K;1:K]); hold on
title("Abs(Average error) vs iteration");grid on;
loglog(1:K,repmat(abs(Bias_analytic(1:N)),1,K));
% Organizing Legend
leg1=[]; leg3=[];
for ii=1:N
    leg1 =[ leg1; strcat("Abs Avg err r_",int2str(ii))];
    leg3 = [leg3; strcat("Abs An. bias r_",int2str(ii))];
end
legend([leg1; leg3],"Location","Southwest");

subplot(2,2,4)
semilogx(1:K,angle(cumsum(err_n,2)./[1:K;1:K])); hold on
title("Phase(Average error) vs iteration");grid on;
semilogx(1:K,repmat(angle(Bias_analytic(1:N)),1,K));
% Organizing Legend
leg2=[]; leg4=[];
for ii=1:N
    leg2 = [leg2; strcat("Phase Avg err r_",int2str(ii))];
    leg4 = [leg4; strcat("Phase An. bias r_",int2str(ii))];
end
legend([leg2; leg4],"Location","Southwest");

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));