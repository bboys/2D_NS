function [setup] = setupAMG(setup, setupType)
if nargin < 1
	setup = [];
end
if nargin < 2
	setup.amg.levels = default('Number of levels for V-cycle AMG', 4);
	setup.amg.maxIt = default('Maximum number of V-cycles', 2);
	setup.amg.tol = default('Relative tolerance', 1e-8);

	setup.amg.nrPreSmooth = default('Number of pre-smoothing steps', 1);
	setup.amg.nrPostSmooth = default('Number of post-smoothing steps', 1);
	setup.amg.smoothType = default({'Type of smoother', 'GS', 'SOR'},1);

	setup.amg.coarseMethod = default({'Type of coarsener', 'Ruge-Stuben', 'PMIS', 'CJLP'},2);
	setup.amg.interpMethod = default({'Type of interpolation', 'Classical', 'F-F'},2);
	setup.amg.theta = default('Connectivity threshold, theta', 0.7);
elseif strcmp(setupType, 'default')
	setup.amg.levels = 4;
	setup.amg.maxIt = 15;
	setup.amg.tol = 1e-8;

	setup.amg.nrPreSmooth = 1;
	setup.amg.nrPostSmooth = 1;
	setup.amg.smoothType = 1;

	setup.amg.coarseMethod = 1;
	setup.amg.interpMethod = 1;
	setup.amg.theta = 0.7;
elseif strcmp(setupType, 'defaultPrecon')
	setup.amg.levels = 4;
	setup.amg.maxIt = 2;
	setup.amg.tol = 1e-8;
	
	setup.amg.nrPreSmooth = 1;
	setup.amg.nrPostSmooth = 1;
	setup.amg.smoothType = 1;

	setup.amg.coarseMethod = 2;
	setup.amg.interpMethod = 2;
	setup.amg.theta = 0.7;
end
end