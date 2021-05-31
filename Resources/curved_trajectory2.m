function [trajectory] = curved_trajectory2(z1,z2,lambdas)
%CURVED_TRAJECTORY1 Returns a curved trajectory between the two numbers
T=length(lambdas);
trajectory=((1-lambdas)*abs(z1)+lambdas*abs(z2)).*exp(1i*((1-lambdas)*angle(z1)+lambdas*(angle(z2)+2*pi)));

end

