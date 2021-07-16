function [SI_sage]=SAGE_biexp4p(x,te)

%x should be [SI_I SI_II R2s R2] where d = SI_I/SI_II; 
tn=te;
TE=te(end);
% % TE = 0.09;
% % TE = 0.087;

SI_sage = zeros(size(tn));
%Create piece-wise function based on Schmiedeskamp H et al. MRM 2012
%67:378-388
SI_sage(tn<(TE/2)) = x(1).*exp(-tn(tn<(TE/2)).*x(3));
SI_sage(tn>(TE/2)) = x(2)*exp(-TE*(x(3)-x(4))).*exp(-tn(tn>(TE/2))*(2*x(4)-x(3)));