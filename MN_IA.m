function [IA_mat] = MN_IA(M_r, k)
% Modified Nested Index assignment used in VaishMDSQ(R, k)
% R: bpss per channel, M_r = sqrt(2.^R);
% k: diagonals (2k+1) of index assignment matrix used
% IA_mat: Index Assignment matrix
q = floor(M_r./k); % Assuming R1=R2
r = mod(M_r,k);                                                                                                                                  
iseven = ~mod(q,2);
if (iseven)
    scan1_end = q/2 - 1;    
else
    scan1_end = (q-3)/2;    
end

IA_mat = [];
for p = 0:scan1_end
for i = 2*p*k+1:(2*p+1)*k;
    ki = min(k,M_r - i);
    IA_mat = [IA_mat;EastScan(i, ki)];
end
for i = (2*p+1)*k+1:(2*p+2)*k;
    ki = min(k,M_r - i);
    IA_mat = [IA_mat;SouthScan(i, ki)];
    % for a 3D IA_mat:
%     IA_mat = cat(3,IA_mat,[SouthScan(i, ki);zeros(2*(k-ki),2)]); 

end
end

if(iseven)
    for(i = M_r-r+1:M_r)
        ki = min(k,M_r - i);
        IA_mat = [IA_mat;EastScan(i, ki)];%         
    end
else
    for(i = (q - 1)*k+1:q*k)
        ki = min(k,M_r - i);
        IA_mat = [IA_mat;EastScan(i, ki)];%         
    end
     for(i = M_r-r+1:M_r)
        ki = min(k,M_r - i);
        IA_mat = [IA_mat;SouthScan(i, ki)];%         
     end
end


% Hafsa Qureshi, www.TSP.ECE.McGill.CA
% 
% $Id: MN_IA.m 2009/06/18 MN_IA-v1.0 $

return

function [E_ki] = EastScan(i,ki)
% index sequence produced during modified east scan, used in MN_IA(M_r, k)
% i: 
% ki: 
% E_ki:

e_1 = [0:ki;zeros(1,ki+1)];
e_1 = e_1(:);
e_v1 = e_1(1:end-1);
e_v2 = e_1(2:end);
E_ki = repmat([i,i],length(e_v1),1)+[e_v1,e_v2];

return

function [S_ki] = SouthScan(i,ki)
% index sequence produced during modified south scan, used in MN_IA(M_r, k)
% i: 
% ki: 
% S_ki:
e_1 = [0:ki;zeros(1,ki+1)];
e_1 = e_1(:);
s_v1 = e_1(2:end);
s_v2 = e_1(1:end-1);
S_ki = repmat([i,i],length(s_v1),1)+[s_v1,s_v2];

return
