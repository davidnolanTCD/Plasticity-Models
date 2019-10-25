function [out]=Celastic(E,nu)
% del = eye(3,3);
% E = 1000; nu = 0.3;
lamda = (nu*E)/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

% for i = 1:3
%     for j = 1:3
%         for k = 1:3
%             for l = 1:3
%                 out(i,j,k,l) = lamda*del(i,j)*del(k,l) + mu*(del(i,k)*del(j,l)...
%                     + del(i,l)*del(j,k));
%             end
%         end
%     end
% end

a = 2*mu + lamda;
out = [ a lamda lamda 0 0 0
    lamda a lamda 0 0 0
    lamda lamda a 0 0 0
    0 0 0 mu 0 0
    0 0 0 0 mu 0
    0 0 0 0 0 mu];


end




