clear all

N = 2^6;
A = gallery('poisson', N);
b = ones(N^2,1);

[setup] = setupAMG([], 'default');
setup.amg.maxIt = 30;
setup.amg.coarseMethod = 1; % PMIS

tic
[amgSystem, setup] = createAmgSystem(A, setup);
fprintf('Setting up amgSystem done in %d seconds\n',toc)
u = zeros(N^2,1);

tic
[u, relres] = amgSolve(A, b, u, setup, amgSystem);
[amgSystem, setup] = createAmgSystem(A, setup);
fprintf('Solving system done in %d seconds\n',toc)

figure
subplot(1,2,1)
semilogy(relres)
subplot(1,2,2)
plot(relres(2:end)./relres(1:end-1))