function [x,hist]=gmresPrecon(A,b,x,setup,precon)
	maxMV = setup.linsolve.maxIt; % max nr of MVs
	tol = setup.linsolve.tol;

    restart = setup.linsolve.gmresRestart;
    restart = min(restart, maxMV);

    kmax = ceil(maxMV/restart); % max nr of outer iterations

	FOM= 0;

    stopCrit= 0; % equals one if converged

	%%%
	[n,m]=size(b); 
    r=b-A*x; 
    nrMV = 1;


	%%% scale tol for relative residual reduction
    rho0=norm(r); rho=rho0; hist=rho;
	tol=rho*tol;

	V=zeros(n,restart); H=zeros(restart+1,restart); gamma=zeros(1,restart);

    for outerk = 1:kmax
    	V(:,1)=r/rho; gamma(1,1)=1; nu=1;
    	for k=1:restart
            if rho<tol, stopCrit= 1; break, end % converged
            [V,H]=ArnStep(A,V,H,k,1, setup, precon);
            nrMV = nrMV + 1;
            %%% compute the nul vector of H'
            gk=(gamma(1,1:k)*H(1:k,k))/H(k+1,k); gamma(1,k+1)=-gk;
            %%% Compute the residual norm
            if FOM
                rho=rho0/abs(gk);
            else
                nu=nu+gk'*gk; rho=rho0/sqrt(nu);
            end
            hist=[hist,rho];
            if nrMV >= maxMV, stopCrit= 1; break, end
    	end % end inner loop
        %%% solve the low dimensional system
        k=size(H,2); e=zeros(k,1); e(1,1)=1;
        if FOM
         y=H(1:k,1:k)\e;
        else
         e=[e;0]; y=H(1:k+1,1:k)\e;
        end

        %%% compute the approximate solution 
        t=  V(:,1:k)*(rho0*y);
        x = x + preconFunc(t, setup, precon);
        % check convergence in inner loop
        if stopCrit| nrMV == maxMV - 1
            break
        else
            % restart
            % r = r - rho0*V(:,1:k)*(H(1:k+1,1:k)*y); % other way of computing b-A*x
            r=b-A*x;  
            nrMV = nrMV + 1;
            rho0 = norm(r); rho=rho0; hist=[hist,rho];
            % tol=rho*tol;
        end

    end % end outer loop


end

function [V,H]=ArnStep(A,V,H,k,kappa, setup, precon)

   if nargin<4, kappa=0.2; end

   w=A*preconFunc(V(:,k), setup, precon);
   [v,h]=Orth(V(:,1:k),w,kappa);
   V(:,k+1)=v;
   H(1:k+1,k)=h;

end

function [v,h]=Orth(V,w,kappa)
%  w=[V,v]*h; with v such that v'*V=0;
%
%  kappa=0: classical Gram-Schmidt
%  kappa>0: repeated Gram-Schmidt with DGKS stopping criterion
%           repeat if tan(angle)<kappa
%  kappa<0: modified Gram-Schmidt
%

global REPEATED

   if nargin==2, kappa=0; end
   [n,k]=size(V);
   if (k>=n), kappa=0; end
   if (k==0), h=norm(w,2); v=w/h; REPEATED=[]; return, end
   if ~exist('REPEATED'), REPEATED=[]; end

   zero=k*k*n*norm(w)*eps; %%% numerical zero 
   v=w;
   if kappa<0 %%% modified Gram-Schmidt
      for j=1:k
        h(j,1)=V(:,j)'*v;
        v=v-V(:,j)*h(j,1);
      end
      rho=norm(v,2);
   else %%% repeated Gram-Schmidt
      h=V'*v;
      v=v-V*h;
      rho=norm(v,2);
      mu=norm(h,2);        t=0;
      %%% Daniel Gragg Kaufman Stewart update criterion
      %%% if kappa==0, classical Gram-Schmidt
      while (rho<kappa*mu & rho>zero)
         g=V'*v; v=v-V*g;
         rho=norm(v,2); mu=norm(g,2);
         h=h+g;           t=t+1;
      end,                REPEATED(1,size(V,2))=t;
   end
   if rho<2*zero
      %%% if full dimension then no expansion
      if (k>=n), v=zeros(n,0); return, end
      %%% if W in span(V) expand with random vector
      v=Orth(V,rand(n,1),kappa); 
      h=[h;0]; 
   else
      h=[h;rho];
      v=v/rho;
   end
end
