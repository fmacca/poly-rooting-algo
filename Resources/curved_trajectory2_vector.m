function [a_n] = curved_trajectory2_vector(a1,a2,lambdas)
%CURVED_TRAJECTORY1_VECTOR Summary of this function goes here
%   Detailed explanation goes here
if length(a1)~=length(a2)
    print("a1 and a2 hve different lenghts");
else
    M=length(a1);
    T=length(lambdas);
    a_n=zeros(M,1,T);
    a_n(1,1,:)=ones(1,1,T);
    for m=2:M
        a_n(m,1,:)=curved_trajectory2(a1(m),a2(m),lambdas);
    end
end

end

