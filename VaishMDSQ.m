function [IA_mat,A,x0_rec,x1_rec,x2_rec,dist] = VaishMDSQ(inp,R,k,md,lamd1,lamd2,del,A,IA_mat)
% [IA_mat,A,x0_rec,x1_rec,x2_rec,dist] = VaishMDSQ(inp,R,k,md,lamd1,lamd2,del,A,IA_mat)
% A Multiple Description Scalar Quantizer based on the paper of V.A.
% Vaishampayan:"Design of Multiple Description Scalar Quantizers",1993
% input: inp,R,k,md,lamd1,lamd2,del,A,IA_mat
% mode: Modified Linear or Modified Nested index assignment
% output: IA_mat,A,x0_rec,x1_rec,x2_rec,dist

% Hafsa Qureshi, www.TSP.ECE.McGill.CA
% $Id: VaishMDSQ.m 2009/06/17 VaishMDSQ-v1.0 $

M = 2.^(2*R);
M_r = sqrt(M); % M_r = M1 = M2 = 2.^R

if nargin~=9
    IA_mat = [];
    if md ==1    
        IA_mat = MN_IA(M_r, k);
    else        
        IA_mat = ML_IA(M_r, k);
    end    
end

N = length(IA_mat(:,1));  % or N can be expressed in terms of M_r & ki

% initial central partition as N cells of equal length (N+1 uniformly
% distributed points) b/w [xmax,xmin] if A_init not specified:
if nargin==7    
    xmax = 1.5;  %instead of max(inp)
    xmin = -1.5; %instead of min(inp)
    A = zeros(N,2);
    A_init = xmin : (xmax-xmin)/(N) : xmax; 
    A(2:end,1) = A_init(2:end-1);
    A(1:end-1,2) = A_init(2:end-1);
    A(1,1) = -inf;
    A(end,2) = inf;    
%     A(1:end,1) = A_init(1:end-1);
%     A(1:end,2) = A_init(2:end);
end

ki = min(k,M_r - [1:M_r]);
l = 1;
Lagr_tol = del+1;
Lagr = inf;
% D1 = 0.2;
% D2 = 0.2;

while (Lagr_tol > del)     
        
    % Determine optimum decoder (21)-(23)
    [x0_rec,x1_rec,x2_rec] = Vaish_dec(IA_mat,A,inp,N,M_r);
    
    % update central partition using extreme point algorithm
    alpha = x0_rec + lamd1*x1_rec([IA_mat(:,1)]) + lamd2*x2_rec([IA_mat(:,2)]);
    beta = x0_rec.^2 + lamd1*x1_rec([IA_mat(:,1)]).^2 + lamd2*x2_rec([IA_mat(:,2)]).^2;
    [al_s,al_si] = sort(alpha);
    be_s = beta(al_si);
    IA_mat(al_si,:) = IA_mat;  % modify Index Assignment mat
    A = CentPart(al_s,be_s,al_si,N); 
    
    % modifying IA matrix
    z_i = find((A(:,1)==0)&(A(:,2)==0));     
    A(z_i,:) = [];
    IA_mat(z_i,:) = [];
    x0_rec(z_i) = [];
    N = length(IA_mat(:,1));
      
    % Compute Lagrangian using (6)
    [d_0,d_1,d_2] = ms_dist(inp,x0_rec,x1_rec,x2_rec,A,IA_mat);    
%     Lagr = [Lagr, d_0 + lamd1*(d_1-D1) + lamd2*(d_2-D2)]; %incorporating D1,D2
    Lagr = [Lagr, d_0 + lamd1*(d_1) + lamd2*(d_2)]; %not incorporating D1,D2
    Lagr_tol = abs(Lagr(end)-Lagr(end-1))./abs(Lagr(end));    
    fprintf(1, 'Cycle %d :   D0= %8.7f,D1= %8.7f,D2= %8.7f\n', l, d_0, d_1, d_2);
    l = l+1;
end
 dist = [d_0,d_1,d_2];
return




 
    
    
    