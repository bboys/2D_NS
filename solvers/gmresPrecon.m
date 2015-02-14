function [x,hist,t]=gmresPrecon(A,b,x,setup,precon,Q,typeSolve)

  if typeSolve == 1 % stokes
    preconType = setup.linsolve.precon;

  elseif typeSolve == 2 % NS
    preconType = setup.nonlin.precon;

  end
  
  kmax = setup.linsolve.maxIt;
  tol = setup.linsolve.tol;
 
  FOM=strcmp(lower(flag),'fom');

  t0=clock; 
  %%%
  [n,m]=size(b); 

  r=b-A*preconFunc(x, setup, precon, Q,preconType);
  rho0=norm(r); rho=rho0; hist=rho;
  %%% scale tol for relative residual reduction
  tol=rho*tol;
  
  V=zeros(n,kmax); H=zeros(kmax+1,kmax); gamma=zeros(1,kmax);
  V(:,1)=r/rho; gamma(1,1)=1; nu=1;
  for k=1:kmax
     if rho<tol, break, end
     [V,H]=ArnStep(A,V,H,k,-1, setup, precon, Q,preconType);
     %%% compute the nul vector of H'
     gk=(gamma(1,1:k)*H(1:k,k))/H(k+1,k); gamma(1,k+1)=-gk;
     %%% Compute the residual norm
     if FOM
        rho=rho0/abs(gk);
     else
        nu=nu+gk'*gk; rho=rho0/sqrt(nu);
     end
     hist=[hist,rho];
  end
  %%% solve the low dimensional system
  k=size(H,2); e=zeros(k,1); e(1,1)=1;
  if FOM
     y=H(1:k,1:k)\e;
  else
     e=[e;0]; y=H(1:k+1,1:k)\e;
  end
  %%% compute the approximate solution 
  x=V(:,1:k)*(rho0*y);

  x = preconFunc(x, setup, precon, Q,preconType);

  t=etime(clock,t0);

return

function [V,H]=ArnStep(A,V,H,k,kappa, setup, precon, Q,preconType)

   if nargin<4, kappa=0.2; end

   w=A*preconFunc(V(:,k), setup, precon, Q,preconType);
   [v,h]=Orth(V(:,1:k),w,kappa);
   V(:,k+1)=v;
   H(1:k+1,k)=h;

return

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
return
