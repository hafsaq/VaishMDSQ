function [IA_mat] = ML_IA(M_r, k)
% Modified Linear Index assignment used in VaishMDSQ(R, k)
% R: bpss per channel, M_r = sqrt(2.^2*R);
% k: diagonals (2k+1) of index assignment matrix used
% IA_mat: Index Assignment matrix
q = floor((M_r-1)/(2*k)); % Assuming R1=R2
r = mod((M_r-1),(2*k));                                                                                                                                  
evenq = ~mod(q,2);
a = 2 + r;
b = M_r + 2*k + 2;
c = M_r + 2;


IA_mat = [];
if (evenq) 
    evenj = 1;
    for j = 2:r+1
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;DownScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;DownScanOj(i, kj)];            
        end 
        evenj = ~evenj;
    end
    for p = 0:q/2-1
        for j = 4*k*p + a : 2*k*(2*p+1)-1+a               
            if (evenj)                
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;UpScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;UpScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
        for j = 2*k*(2*p+1)+a : 2*k*(2*p+2)-1+a   
            if (evenj)
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;DownScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;DownScanOj(i, kj)];
            end 
            evenj = ~evenj;
        end
    end   
    IA_mat = [IA_mat;[M_r/2, M_r/2+1]];   
    evenj = 1;
    for p = 0:q/2-1
        for j = 4*k*p + c : 2*k*(2*p+1)-1+c            
            if (evenj)                
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;UpScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;UpScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
        for j = 2*k*(2*p+1)+c : 2*k*(2*p+2)-1+c              
            if (evenj)
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;DownScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;DownScanOj(i, kj)];
            end 
            evenj = ~evenj;           
        end
    end 
    evenj = 1;
    for j = 2*M_r-r+1 : 2*M_r
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;UpScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;UpScanOj(i, kj)];
        end                    
        evenj = ~evenj;
    end
    
  %  evenj toggle, scans done till this part
else % q odd
    evenj = 1;
    for j = 2:r+1
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;UpScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;UpScanOj(i, kj)];
        end                    
        evenj = ~evenj;
    end
    for p = 0:(q-3)/2
        for j = 4*k*p + a : 2*k*(2*p+1)-1+a 
            if (evenj)                 
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;DownScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;DownScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
        for j = 2*k*(2*p+1)+a : 2*k*(2*p+2)-1+a
            if (evenj)                
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;UpScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;UpScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
    end   
    evenj = ~mod(M_r-2*k+1,2);
    for j = M_r-2*k+1:M_r
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;DownScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;DownScanOj(i, kj)];
        end                    
        evenj = ~evenj;
    end
    IA_mat = [IA_mat;[M_r/2, M_r/2+1]];    
    evenj = 1;
    for j = M_r+2:M_r+1+2*k
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;UpScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;UpScanOj(i, kj)];
        end                    
        evenj = ~evenj;
    end
    for p = 0:(q-3)/2
        for j = 4*k*p + b : 2*k*(2*p+1)-1+b 
            if (evenj)                 
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;DownScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;DownScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
        for j = 2*k*(2*p+1)+b : 2*k*(2*p+2)-1+b
            if (evenj)                 
                i = j/2;
                kj = min([i-1,floor(k/2),M_r-i]);
                IA_mat = [IA_mat;UpScanEj(i, kj)];
            else
                i = (j-1)/2;
                kj = min([i-1,floor((k-1)/2),M_r-i-1]);
                IA_mat = [IA_mat;UpScanOj(i, kj)];
            end             
            evenj = ~evenj;
        end
    end
    evenj = 1;
    for j = 2*M_r-r+1:2*M_r
        if (evenj)
            i = j/2;
            kj = min([i-1,floor(k/2),M_r-i]);
            IA_mat = [IA_mat;DownScanEj(i, kj)];
        else
            i = (j-1)/2;
            kj = min([i-1,floor((k-1)/2),M_r-i-1]);
            IA_mat = [IA_mat;DownScanOj(i, kj)];
        end                    
        evenj = ~evenj;
    end
end

% Hafsa Qureshi, www.TSP.ECE.McGill.CA
% 
% $Id: MN_IL.m 2009/08/08 ML_IA-v1.0 $

return

function [U_kj] = UpScanEj(i,kj)
% index sequence produced during modified up scan, used in ML_IA(M_r, k)
% i: 
% kj: 
% U_kj:

v2 = -kj:kj;
v1 = fliplr(v2);
U_kj = [v1.',v2.'] + i;
return
function [U_kj] = UpScanOj(i,kj)
% index sequence produced during modified up scan, used in ML_IA(M_r, k)
% i: 
% kj: 
% U_kj:

v2 = -kj:kj+1;
v1 = fliplr(v2);
U_kj = [v1.',v2.'] + i;
return
function [D_kj] = DownScanEj(i,kj)
% index sequence produced during modified down scan, used in ML_IA(M_r, k)
% i: 
% kj: 
% D_kj:

v1 = -kj:kj;
v2 = fliplr(v1);
D_kj = [v1.',v2.'] + i;
return
function [D_kj] = DownScanOj(i,kj)
% index sequence produced during modified down scan, used in ML_IA(M_r, k)
% i: 
% kj: 
% D_kj:

v1 = -kj:kj+1;
v2 = fliplr(v1);
D_kj = [v1.',v2.'] + i;
return