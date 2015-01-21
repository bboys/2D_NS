load testSolveBig.mat
b = rhsVec;

tic
x = A\b;
toc
norm(b - A*x)/norm(b)





tic
% setup.type='nofill'; 
% [L,U] = ilu(A,setup);
L = []; U = [];
[~, flag, relres, resvec] = gmres(A, b, [], 10^-10, 500, L, U)
toc

