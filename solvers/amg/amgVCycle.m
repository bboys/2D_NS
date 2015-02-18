function [uh, res, amgSystem] = amgVCycle(uh, fh, level, maxLvl, amgSystem, setup)
% solve amgSystem.level(1)*u = f
% maxLvl means, maxLvl = 1 <=> 2 cycle AMG
% current level: start with 1, end at maxLvl

% pre smooth
[uh, amgSystem] = amgSmoother(uh, level, fh, setup, amgSystem, setup.amg.nrPreSmooth);


rh = fh - amgSystem.level(level).matrix*uh;

% always one step behind (but does not cost MV)
if level == 1,
	res = norm(rh);
end

% restrict residual
r2h = amgSystem.level(level).interp'*rh;

% recursive call or solve at coarsest level = maxLvl
if level == maxLvl
	% e2h = amgSolve(r2h, level, amgSystem, setup)
	e2h = amgSystem.level(level + 1).matrix\r2h;
else
	e2h = amgVCycle(zeros(size(r2h)), r2h, level + 1, maxLvl, amgSystem, setup);
end

% interpolate error
eh = amgSystem.level(level).interp*e2h;
	
% update solution
uh = uh + eh;

% post smooth
[uh, amgSystem] = amgSmoother(uh, level, fh, setup, amgSystem, setup.amg.nrPostSmooth);

end