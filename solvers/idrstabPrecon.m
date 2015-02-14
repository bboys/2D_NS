function [x,resvec,flag] = idrstabPrecon(A,b,s,l,Rs,x,setup,precon,Q)
% Implementation of the IDR(s)stab(l) algorithm, as proposed by G. Sleijpen
% and M van Gijzen in 'Exploiting BiCGstab(l) strategies to induce
% dimension reduction' (2010).
%
% Input/output given as
%
% [x,resvec,flag] = IDRstab(A,b,tol,maxit,s,l,Rs,M1,M2,x)
%
% Reduce the residual gap by replacing all primary residuals,
% as proposed by K. Aihara (2014).
%
% Rs is n x s shadow residual
% Rm = (rm,Arm,...) n x (l+1)
% R = (r,Ar,...) n x (l+1)
% U = (Um,AUm,...) n x (s*(l+1))
% V = (V,AU,...) n x (s*(l+1))

% Check input variables
n = length(b);

% if nargin < 10 || isempty(x), x = zeros(n,1); end

% if nargin < 5 || isempty(s), s = 4; end
if isempty(Rs), Rs = orth(randn(n,s)); end
% if nargin < 6 || isempty(l), l = 2; end
% if nargin < 4 || isempty(maxit), maxit = 2*n; end
% if nargin < 3 || isempty(tol), tol = 10^-8; end
% if nargin < 2, display('Not enough input arguments'), end

tol = setup.linsolve.tol;
maxit = setup.linsolve.maxIt;


% warning off % Suppress mldivide warnings

% Initialization
flag = 0;
n = length(b);
R = zeros(n,l+1);
R(:,1) = b -  A*preconFunc(x, setup, precon, Q);
U = zeros(n,s*(l+1));
V = zeros(n,s*(l+1));
iter = 1;
resvec = zeros(maxit+1,1);
resvec(1) = norm(R(:,1));
tol = tol*norm(b);
sigma = eye(s); 

% Precompute Rs'*A
Rsold = Rs';
Rs = Rsold*A;
Rs = preconFunc(Rs', setup, precon, Q)';

resvec(iter+1:s+iter) = resvec(1)*ones(s,1);
iter = iter + s;


while resvec(iter)>tol
    % The IDR step
    for j = 1:l
        % Iteratively build basis and residual inside G_{k+j}'
        if j == 1 
            % If rrep is used then for the first projection we need the old
            % shadow residual
            alpha = sigma\(Rsold*R(:,1));
        else
            alpha = sigma\(Rs*R(:,j-1)); 
        end
        v1 = U(:,1:s)*alpha;
        x = x + v1;
        Av =  A*preconFunc(v1, setup, precon, Q);

        % If rgap then residual is replaced (using MV) instead of
        % projected
        R(:,1) = R(:,1) - Av;
        iter = iter + 1;
        resvec(iter) = norm(R(:,1));  

        for i = 1:j-2
            R(:,i+1) = R(:,i+1) - U(:,(i+1)*s+1:(i+2)*s)*alpha;
        end

        if j > 1
            % Building A^(j-1)*r
            v = R(:,j-1);

            R(:,j) = A*preconFunc(v, setup, precon, Q);

            iter = iter + 1;
            resvec(iter) = norm(R(:,1));
        end

        % Build basis for G_{k+j}': A^i*V*eq, for q = 1:s, i =0:j+1-rgap
        for q = 1:s
            % Projections into next space...
            if q == 1
                rr = Rs*R(:,j);
                alpha = sigma\rr;
                for i = 0:j-1
                    V(:,i*s+1) = R(:,i+1) - U(:,(i)*s+1:(i+1)*s)*alpha;        
                end
            else 
                rr = Rs*V(:,(j)*s+q-1);
                alpha = sigma\rr;
                for i = 0:j-1
                    V(:,i*s+q) = V(:,(i+1)*s+q-1) - U(:,i*s+1:(i+1)*s)*alpha;
                end
            end

            % Expanding
            v = V(:,(j-1)*s+q);
 
            V(:,(j)*s+q) = A*preconFunc(v, setup, precon, Q);

            % Orthogonal basis
            for k = 1:q-1
                mu = (V(:,(j)*s+k)'*V(:,(j)*s+q))/(V(:,(j)*s+k)'*V(:,(j)*s+k));
                for i = 0:j
                    V(:,i*s+q) = V(:,i*s+q) - V(:,i*s+k)*mu;
                end 
  
            end
            % Normalisation
            cnorm = norm(V(:,(j)*s+q));
            for i = 0:j
                V(:,i*s+q) = V(:,i*s+q)/cnorm;
            end  
             
            iter = iter + 1;
            resvec(iter) = norm(R(:,1));

            if iter > maxit || resvec(iter) < tol,break,end    
        end
        if iter > maxit || resvec(iter) < tol,break,end
        
        if j == l
            % If j = l then we need A^l*r in order to perform the
            % polynomial step
            v = R(:,j);

            v = preconFunc(v, setup, precon, Q);
            v = A*v;

            R(:,j+1) = Av;    
            iter = iter + 1;
            resvec(iter) = norm(R(:,1));  
        end
        U = V;
        if j < l
            sigma = Rs*V(:,(j)*s+1:(j+1)*s);
        end
    end
    if iter > maxit ,flag = 1;break, elseif resvec(iter) < tol,break,end

    % Polynomial step: minize residual (enter G_{k+l})
    gamma = R(:,2:l+1)\R(:,1);
    v = R(:,1:l)*gamma;

    % Reducing residual gap...
    Av =  A*preconFunc(v, setup, precon, Q);
    R(:,1) = R(:,1) - Av;   
    iter = iter + 1;
    resvec(iter) = norm(R(:,1)); 
    x = x + v;

    % Prepare for next IDR step
    for j = 1:l
        U(:,1:s) = U(:,1:s) - gamma(j)*U(:,(j)*s+1:(j+1)*s);
    end  
    sigma = Rs*U(:,1:s);
end
resvec = resvec(1:iter);  
x = preconFunc(x, setup, precon, Q);
end

