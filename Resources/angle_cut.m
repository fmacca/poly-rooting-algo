function [arg] = angle_cut(z,cut)
%ANGLE_CUT Returns the argument of z, the angle retuned is in the interval
%[-2*pi+cut, cut].
%cut must be in [-pi, pi]

arg=angle(z);
arg(arg>cut)=arg(arg>cut)-2*pi;


end

