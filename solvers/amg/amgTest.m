
profile on

N = 16;
A = gallery('poisson', N);
b = ones(N^2,1);

[setup] = setupAMG([], 'default');
setup.amg.maxIt = 30;
setup.amg.coarseMethod = 2; % RS, PMIS
setup.amg.interpMethod = 2; % classical, F-F
setup.amg.smoothType = 1; % GS, symmetric GS
setup.amg.levels = 7;
setup.amg.theta = 0.7;

tic
[amgSystem, setup] = createAmgSystem(A, setup);
fprintf('Setting up amgSystem done in %d seconds\n',toc)
u = zeros(N^2,1);

tic
[u, relres] = amgSolve(A, b, u, setup, amgSystem);
fprintf('Solving system done in %d seconds\n',toc)

profile viewer

figure
subplot(1,2,1)
semilogy(relres)
subplot(1,2,2)
plot(relres(2:end)./relres(1:end-1))

% figure
% for level = 1:setup.amg.levels + 1
% 	subplot(2,3,level)
% 	spy(amgSystem.level(level).matrix)
% 	title(sprintf('Theta = %3.2f, fill = %4.2f %%', [setup.amg.theta,...
% 		100*nnz(amgSystem.level(level).matrix)/...
% 		prod(size(amgSystem.level(level).matrix))]))
% end