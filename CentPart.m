function A = CentPart(alpha,beta,al_si,N)
% Cantral Partitioning, used in VaishMDSQ(arg1, arg2)
% N: cardinality of central partition set
% alpha: sorted alpha
% beta: beta corresponding to sorted alpha
% A : lower and upper end-points of the central partitions in 1st & 2nd
% columns resp.

m = 1;
p = m + 1;
t_L = zeros(1,N);
t_U = zeros(1,N);
t_L(m) = -inf;
n = N;
h_l = [];

while (1)
    t = (beta(m)-beta(p))./(2*(alpha(m)-alpha(p)));
    l = m+1 : n;   
    h_l = 2*alpha(l).*t - beta(l);
    p_bar = find([h_l==max(h_l)]) + m;
    if length(p_bar)>1
        p_bar = p_bar(end);
    end
    if(p_bar~=p)
        p = p_bar;
        n = p_bar;
    else
        t_L(p) = t;
        t_U(m) = t;
        if (p == N)
            break;
        else
            m = p;
            p = m + 1;
            n = N;
        end
    end
    h_l = [];
end
t_U(N) = inf;
% t_L_mod(al_si) = t_L;
% t_U_mod(al_si) = t_U;
% A = [t_L_mod.',t_U_mod.'];
A = [t_L.',t_U.'];

% Hafsa Qureshi, www.TSP.ECE.McGill.CA
% 
% $Id: VaishMDSQ.m 2009/06/18 VaishMDSQ-v1.0 $

return