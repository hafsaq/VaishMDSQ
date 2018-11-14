function [d_0,d_1,d_2] = ms_dist(inp,x0_rec,x1_rec,x2_rec,A,IA_mat)
% DESCRIPTION:
%       
% INPUT:
%       
% OUTPUT: 
%       
% USAGE:
%       
% AUTHOR:
%       Hafsa Qureshi
% HISTORY: 
%     08/04/09:  Created


ind = [];
inp_l = length(inp);
[n,bin] = histc(inp,[A(:,1);A(end,2)]);
d_0=0; d_1=0; d_2=0; 


for i=1:length(x0_rec)
    ind = find(bin==i);
    d_0 = d_0 + sum((inp(ind) - x0_rec(i)).^2)./inp_l;    
    d_1 = d_1 + sum((inp(ind) -  x1_rec(IA_mat(i,1))).^2)./inp_l;
    d_2 = d_2 + sum((inp(ind) -  x2_rec(IA_mat(i,2))).^2)./inp_l;
end


return