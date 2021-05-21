clear all
close all
clc
%% Parameters
N=4; % order of the polynomial
sigma_a=.001; % variance of the noise on coefficients
n_simulations=20; % number of simulations: n of times we choose a polynomial and a coavariance matrix
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
plot_layout="big"; %"side" % a parameter setting layout for the plot can be "big" or "side"

%% Useful variables
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Matrix for change of variables

%% Simulation
for simul = 1:n_simulations
    disp(strcat('Polynomial n.',num2str(simul)));
    % Generate circular/proper covariance matrix
    flag=0; %flag to become 1 when a valid (spd) matrix is generated
    while(~flag)
        Temp = sigma_a*1/2/N*wishrnd(eye(2*N),2*N);
        Sigma=J*Temp*J';
        Gamma=Sigma(1:N,1:N);
        C_atilda=1/4*J'*[Gamma zeros(N);zeros(N) conj(Gamma)]*J;
        try A=chol(C_atilda)';
            disp('Matrix is symmetric positive definite.')
            flag=1;
        catch ME
            disp('Matrix is not symmetric positive definite')
        end
    end
    Sigma=J*C_atilda*J';
    clear Temp Gamma
    % Generate random polynomial
    r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)]; % Generate random roots
    a=conj(poly(r)'); % corresponding noise-free polynomial
    deriv=polyder(a); % Coefficients of polinomial derivarive
    % Simulation
    h = waitbar(0,'Simulations in progress ... Please wait...');
    avgs=zeros(N,1);
    for k=1:K
        noise_tilda=A*randn(2*N,1);
        a_n(:,k)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)];   % coff. a_0=1
        r_n(:,k)=roots(a_n(:,k));
        r_n_ord(:,k)=r_n(Order_roots(r_n(:,k),r),k); % Order the estimated roots
        err_n(:,k)=r_n_ord(:,k)-r; % save the estimation error of roots
        
        waitbar(k/K)
    end
    avgs = mean(r_n_ord,2); % mean of estimated roots
    close(h);
    % Plot of the coefficients distributions
    figure(1);
    if plot_layout=="big"
        subplot(1,2,1);
    else
        subplot(2,1,1);
    end
    [hz,hp,ht] =zplane(1); hold on; % Unitary circle
    for ii=2:N+1
        plot(real(a_n(ii,:)),imag(a_n(ii,:)),'.','MarkerSize',1); hold on; % Simulated coeff
    end
    for ii=1:N
        Ellipse_plot(0.1*inv(C_atilda([ii N+ii],[ii N+ii])),[real(a(1+ii)),imag(a(1+ii))])
    end
    plot(real(a),imag(a),'*k','MarkerSize',20);
    axis(5*[-1,1,-1,1]);
    title("Coefficients");
    % Plot the roots distributuions
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
    plot(real(avgs),imag(avgs),'.b','MarkerSize',15); % Mean estimate
    plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
    axis((scale+1)*[-1,1,-1,1]);
    title("Roots");
    % Plot the error vs iterations
    figure(2);
    subplot(2,1,1)
    loglog(1:K,abs(cumsum(err_n,2))./repmat(1:K,N,1));
    hold on
    subplot(2,1,2)
    loglog(1:K,cumsum(abs(err_n),2)./repmat(1:K,N,1));
    hold on
    % Clear some variables
    clear h a_n err_n r_n r_n_ord
end

%% Save figures
figure(1);
saveas(gcf,'./Circular bias pictures/Nequal4/Gauss plane');
saveas(gcf,'./Circular bias pictures/Nequal4/Gauss plane.png');

figure(2);
saveas(gcf,'./Circular bias pictures/Nequal4/Bias vs iter');
saveas(gcf,'./Circular bias pictures/Nequal4/Bias vs iter.png');
