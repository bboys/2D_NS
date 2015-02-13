N = 5;
A = gallery('poisson', N);
b = ones(N^2,1);

levels = 1;
setup.amg.nrPreSmooth = 3;
setup.amg.nrPostSmooth = 3;

setup.amg.coarseMethod = 'RS';
setup.amg.interpMethod = 'classical';
setup.amg.theta = 0.7;

[amgSystem] = createAmgSystem(A, levels, setup);
u = zeros(N^2,1);

for iter = 1:10
	[u] = amgVCycle(u, b, 1, levels, amgSystem, setup);
	norm(b-A*u)/norm(b)
end