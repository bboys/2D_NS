function [uh, res] = amgVCycle(uh, fh, level, levels, amgSystem, setup)
% solve amgSystem.level(1)*u = f
% levels means, levels = 1 <=> 2 cycle AMG
% current level: start with 1, end at levels

% pre smooth
for nu1 = 1:setup.amg.nrPreSmooth
	[uh, amgSystem] = amgSmoother(uh, level, fh, setup, amgSystem);
end

rh = fh - amgSystem.level(level).matrix*uh;

% always one step behind (but does not cost MV)
if level == 1,
	res = norm(rh);
end

% restrict residual
r2h = amgSystem.level(level).interp'*rh;

% recursive call or solve at coarsest level = levels
if level == levels
	% e2h = amgSolve(r2h, level, amgSystem, setup)
	e2h = amgSystem.level(level + 1).matrix\r2h;
else
	e2h = amgVCycle(zeros(size(r2h)), r2h, level + 1, levels, amgSystem, setup);
end

% interpolate error
eh = amgSystem.level(level).interp*e2h;
	
% update solution
uh = uh + eh;

% post smooth
for nu2 = 1:setup.amg.nrPostSmooth
	[uh, amgSystem] = amgSmoother(uh, level, fh, setup, amgSystem);
end

end