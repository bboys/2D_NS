clear all

N = 5;
A = gallery('poisson', N);
b = ones(N^2,1);

levels = 3;
setup.amg.nrPreSmooth = 1;
setup.amg.nrPostSmooth = 1;

setup.amg.coarseMethod = 'RS';
setup.amg.interpMethod = 'classical';
setup.amg.theta = 0.8;

[amgSystem] = createAmgSystem(A, levels, setup);
u = zeros(N^2,1);

maxIt = 10;

for iter = 1:maxIt
	u = amgVCycle(u, b, 1, levels, amgSystem, setup);
	% u = amgSmoother(u, A, b, setup);
	res(iter) = norm(b-A*u)/norm(b);
end
semilogy(res)