clear all
close all
clc
%% Parameters
N=2; % order of the polynomial
sigma_a=.001; % variance of the noise on coefficients
% Generate a random circular covariance
Temp = sigma_a*1/2/N*wishrnd(eye(2*N),2*N);
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Matrix for change of variables
Sigma=J*Temp*J';
Gamma=Sigma(1:N,1:N);
C_atilda=1/4*J'*[Gamma zeros(N);zeros(N) conj(Gamma)]*J;
try A=chol(C_atilda)';
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end
% A=chol(C_atilda)';
plot_layout="big"; %"side

%% Set the "true" polynomial
scale=2; % the roots will be in Re(r),Im(r) in [-scale +scale]^2, a square
r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
a=conj(poly(r)'); % polynomial noise-free

%% Computing some useful quantities
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Matrix for change of variables
deriv=polyder(a); % Coefficients of polinomial derivarive
Sigma=J*C_atilda*J';

%% Simulation
h = waitbar(0,'Simulations in progress ... Please wait...');
K=10^4; % Number of iterations
avgs=zeros(N,1);
for k=1:K
    noise_tilda=A*randn(2*N,1);
    a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)];   % coff. a_0=1
    % a_n=a+[sigma_a*(randn(N+1,1)+1i*randn(N+1,1))]; %  coff. a_0 rnd
    r_n(:,k)=roots(a_n(:,k));
    r_n_ord(:,k)=r_n(order_roots_permutations(r_n(:,k),r),k);
    err_n(:,k)=r_n_ord(:,k)-r;
%     for ii=1:N
%         [~,indx]=min(abs(r_n(:,k)-ones(N,1)*r(ii)));
%         avgs(ii)=avgs(ii)+r_n(indx,k);%this can take a value twice
%         % try ordering with permutations
%     end
    waitbar(k/K)
end
% avgs=avgs/K;
avgs = mean(r_n_ord,2);
close(h);

%% Plot the results
figure(1);
if plot_layout=="big"
    subplot(1,2,1);
else
    subplot(2,1,1);
end
[hz,hp,ht] =zplane(1); hold on; % Unitary circle
%viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.5);hold on;
for ii=2:N+1
    plot(real(a_n(ii,:)),imag(a_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated coeff
end

for ii=1:N
    ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
end
plot(real(a),imag(a),'*k','MarkerSize',20);

axis(5*[-1,1,-1,1]);
title("Coefficients");hold off

if plot_layout=="big"
    subplot(1,2,2);
else
    subplot(2,1,2);
end
[hz,hp,ht] =zplane(1); hold on; % Unitary circle

for ii=1:N
    plot(real(r_n_ord(ii,:)),imag(r_n_ord(ii,:)),'.','MarkerSize',1); hold on; % Simulated roots
end

% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(r(ii)),imag(r(ii))])
% end

plot(real(avgs),imag(avgs),'.b','MarkerSize',15); % Biased mean estimate

% for ii=1:N
%     plot(real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii)),'.g','MarkerSize',15); % Estimate - bias
% end

plot(real(r),imag(r),'*k','MarkerSize',20); % True roots

axis(3*[-1,1,-1,1]);
title("Roots");hold off

% figure(2);
% corrplot([real(a_n(2:end,:)') imag(a_n(2:end,:)')]);
% 
% title("Correlation between the two");hold off

%% Plot the error
% figure();
% plot(1:K,abs(cumsum(err_n,2))./[1:K;1:K]);
% 
% figure();
% plot(1:K,cumsum(abs(err_n),2)./[1:K;1:K]);

%%
figure();
subplot(2,1,1)
loglog(1:K,abs(cumsum(err_n,2))./[1:K;1:K]);

subplot(2,1,2)
loglog(1:K,cumsum(abs(err_n),2)./[1:K;1:K]);

% figure();
% loglog(1:K,cumsum(abs(err_n),2)./[1:K;1:K]);