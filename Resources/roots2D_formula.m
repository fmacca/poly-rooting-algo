function [r] = roots2D_formula(a,cut)
%ROOTS2D_FORMULA This function gets as input a polynomial of degree 2
%It returns its roots comuted by the formula
%The order is provided by the "cut" (TODO)
if length(a)~=3
    print("Polynomial does not have degree 2");
else
    [a_transf,t] = poly2D_affinetransf(a);
    delta=-a_transf(3);
    r=[t+sqrt(abs(delta))*exp(1i*angle_cut(delta,cut)/2); t+sqrt(abs(delta))*exp(1i*angle_cut(delta,cut)/2+1i*pi)];
end

end

