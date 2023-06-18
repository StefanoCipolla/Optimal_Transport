function NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_vars,pivot_thr,Struct)
% ==================================================================================================================== %
% Newton_factorization: Factorize the Newton matrix
% -------------------------------------------------------------------------------------------------------------------- %
% NS = Newton_factorization(A,A_tr,Q,x,z,delta,rho,pos_vars,free_vars) returns a MATLAB struct that holds the
%      factorization of the Newton matrix for solving the step equations in
%      the IPM, as well as relevant information concerning failure.
% Factorization Method
  % --------------------
  % 1: augmented system, LDL' factorization.
% 
% Authors: S.Cipolla, J. Gondzio, F. Zanetti.
% ==================================================================================================================== %
[m, n] = size(A);
NS = struct();
% ==================================================================================================================== %
% LDL' factorization of KKT matrix
% -------------------------------------------------------------------------------------------------------------------- %
%
% MATLAB uses MA57 when K is sparse, which is not available in OCTAVE. 
% -------------------------------------------------------------------------------------------------------------------- %
NS.x = x;
NS.z = z;
Q_bar = zeros(n,1);
if (size(pos_vars,1) > 0)
    Q_bar(pos_vars) = z(pos_vars)./x(pos_vars) + rho;
    Q_bar(free_vars) = rho;
else
    Q_bar(:) = rho;
end

if strcmp(Struct.Fact,'ldl')
    K = [Q+spdiags(Q_bar,0,n,n), A_tr; A, -spdiags(delta.*ones(m,1),0,m,m)];  
    [NS.L,NS.D,NS.pp] = ldl(K,pivot_thr,'vector'); %Small pivots allowed, to avoid 2x2 pivots.

elseif strcmp(Struct.Fact,'pcg')
    NS.A = A;
    NS.D = 1./Q_bar;  
    numbers   = 1:n;
    index     = numbers(NS.D > 0.4*(Struct.mu/(1+rho*Struct.mu)));
    if index(end)<n
        index     = [index n];
    end
    pDiag         = sparse(index,ones(length(index),1),NS.D(index));
    NS.D_a = pDiag;
    K      = A*(pDiag.*A_tr)+spdiags(delta.*ones(m,1),0,m,m);
    NS.K   = K;
    NS.L = ichol(K,struct('type','ict','droptol',1e-3));
    NS.delta = delta;
    NS.rho = rho;

elseif strcmp(Struct.Fact,'thr+direct')
    NS.A = A;
    NS.D = 1./Q_bar;  
    numbers   = 1:n;
    index     = numbers(NS.D > 0.4*(Struct.mu/(1+rho*Struct.mu)));
    if index(end)<n
        index     = [index n];
    end
    pDiag         = sparse(index,ones(length(index),1),NS.D(index));
    NS.D_a = pDiag;
    K      = A*(pDiag.*A_tr)+spdiags(delta.*ones(m,1),0,m,m);
    NS.K   = K;
    [NS.L,~,NS.pp] = chol(K,'upper','vector');
    NS.delta = delta;
    NS.rho = rho;
else  
    NS.A = A;
    NS.D = 1./Q_bar;
    K    = A*(NS.D.*A_tr)+spdiags(delta.*ones(m,1),0,m,m);
    [NS.L,~,NS.pp] = chol(K,'vector');    
 end
% ==================================================================================================================== %  
 
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
end
