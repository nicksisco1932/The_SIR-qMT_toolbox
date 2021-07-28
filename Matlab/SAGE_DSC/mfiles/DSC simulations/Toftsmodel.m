function F = Toftsmodel(beta,x_array)


global AIFfull;
Ktrans = beta(1);
Ve = beta(2);

fval = Ktrans*conv(AIFfull,exp((-Ktrans/Ve)*x_array));
F = fval(1:size(x_array,2));







% for it = 1:1:size(x_array),
%     t = x_array(it);
%     tempsum =  0.0;
%   
%     for it2 = 1:1:(it-1),
%         t2 = x_array(it2);
%         F2 = @(t2)(AIFfull(it2))*exp((-Ktrans/Ve)*(t-t2));
%         tempsum = quad(F2,0,t2);
%     end;
%  
%   fval(it) = Ktrans*tempsum;
%   
% end;
%size(fval)
%F = fval';