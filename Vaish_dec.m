function [x0_rec,x1_rec,x2_rec] = Vaish_dec(IA_mat,A,inp,N,M_r)
% DESCRIPTION:
%       Optimal decoder given by (21)-(23)
% INPUT:
%       
% OUTPUT: 
%       
% USAGE:
%       [x0_rec,x1_rec,x2_rec] = Vaish_dec(IA_mat,A,inp,N,M_r);
% AUTHOR:
%       Hafsa Qureshi
% HISTORY: 
%     07/13/09:  Created
x0_rec = zeros(1,N);
x1_rec = zeros(1,M_r);
n1 = zeros(1,M_r);
x2_rec = zeros(1,M_r);
n2 = zeros(1,M_r);
ind = [];
[n,bin] = histc(inp,[A(:,1);A(end,2)]);

for i = 1:N
   ind = find(bin==i);     
   x0_rec(i) = mean(inp(ind));
   x1_rec(IA_mat(i,1)) = sum([x1_rec(IA_mat(i,1)),inp(ind)]);
   n1(IA_mat(i,1)) = n1(IA_mat(i,1)) + n(i);
   x2_rec(IA_mat(i,2)) = sum([x2_rec(IA_mat(i,2)),inp(ind)]); 
   n2(IA_mat(i,2)) = n2(IA_mat(i,2)) + n(i);
end
x1_rec = x1_rec./n1;
x2_rec = x2_rec./n2;




return