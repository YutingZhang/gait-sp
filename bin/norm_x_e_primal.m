function [x, e, nIter, nIter_in_total] = norm_x_e_primal(y, A, tol, tol_int, maxIter, maxIter_apg)

% This PALM program solves the problem:
% min \|x\|_1 + \|e\|_1    s.t.      y = Ax + e

% Question? --Zihan Zhou (zzhou7@illinois.edu)
% Copyright: Perception and Decision Laboratory, University of Illinois,
% Urbana-Champaign

[m,n] = size(A) ;

% for real face data, choose tol = 1e-1, tol_int = 1e-5 ~ 1e-6 depending on
% the accuracy/speed trade off
% for synthetic random data, choose tol = 1e-6, tol_int = 1e-3 
if ~exist('tol', 'var')
    tol = 5e-2;
end
if ~exist('tol_int', 'var')
    tol_int = 1e-6;
end
if ~exist('maxIter', 'var')
    maxIter = 200;
end
if ~exist('maxIter_apg', 'var')
    maxIter_apg = 400;
end

%tol_int = 1e-6 ;

VERBOSE = 0;

% G = A'*A ;
opts.disp = 0;
tau = eigs(A*A',1,'lm',opts);
tauInv = 1/tau ;

nIter = 0 ;

mu = 2 *m / norm(y,1);

if VERBOSE
    disp(['mu is: ' num2str(mu) ' tau is: ' num2str(tau)]);
end

lambda = zeros(m,1);
x = zeros(n,1) ;
e = y;

converged_main = 0 ;
nIter_in_total = 0 ;

%maxIter = 50 ;
%maxIter = 1;
%maxIter_apg = 50 ;
%maxIter_apg = 3;

while ~converged_main
    
    muInv = 1/mu ;
    lambdaScaled = muInv*lambda ;
    
    nIter = nIter + 1 ;
    
    e_old_main = e ;
    x_old_main = x ;
    
        
    temp2 = y + lambdaScaled ;
    
    temp = temp2 - A*x ;
    
    %e = shrink(temp, muInv) ;
    e = sign(temp) .* max(abs(temp) - muInv, 0);
    
    converged_apg = 0 ;
    
    temp1 = A'*(e - temp2) ;
    
    nIter_apg = 0 ;
    
    t1 = 1 ; z = x ;
    
    muTauInv = muInv*tauInv ;
    
    Gx = A' * (A * x);
    Gz = Gx;

    while ~converged_apg
            
        nIter_apg = nIter_apg + 1 ;
            
        x_old_apg = x ;
        Gx_old = Gx;
            
        temp = z - tauInv*(temp1 + Gz) ;           
            
        %x = shrink(temp, muTauInv) ;
        x = sign(temp) .* max(abs(temp)-muTauInv, 0);
        
        Gx = A' * (A * x);

        s = tau * (z - x) + Gx - Gz;
        if norm(s) < tol_int * tau * max(1,norm(x))
            converged_apg = 1;
        end
            
        if nIter_apg >= maxIter_apg
            converged_apg = 1 ;
        end
            
        t2 = (1+sqrt(1+4*t1*t1))/2 ;
        z = x + ((t1-1)/t2)*(x-x_old_apg) ;
        Gz = Gx + ((t1-1)/t2) * (Gx - Gx_old);
        t1 = t2 ;
        
    end
    
    nIter_in_total = nIter_in_total + nIter_apg;

    
    if VERBOSE
        disp(['Iteration ' num2str(nIter)]) ;
        disp([norm(x_old_main - x)/norm(x_old_main) norm(e_old_main - e)/norm(e_old_main) nIter_apg])
        
        figure(1);
        subplot(3,1,1);
        plot(x);
        title('x');
        subplot(3,1,2);
        plot(e);
        title('e');
        subplot(3,1,3);
        plot(lambda);
        title('lambda');
        pause(.1);
    end 
        
    lambda = lambda + mu*(y - A*x - e) ;   
    
    if norm(x_old_main - x) < tol * norm(x_old_main) && ...
        norm(e_old_main - e) < tol * norm(e_old_main) 
        converged_main = 1 ;
    end
    
    if ~converged_main && nIter >= maxIter
        disp('Maximum Iterations Reached') ;
        converged_main = 1 ;
    end
    
end