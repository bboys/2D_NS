function [u, amgSystem] = amgSmoother(u, level, f, setup, amgSystem, repeat)
if setup.amg.smoothType == 1 % GS

	for rep = 1:repeat
		u = amgSystem.level(level).Ls\(f - amgSystem.level(level).U*u);
	end
elseif setup.amg.smoothType == 2 % symmetric GS

	for rep = 1:repeat
		u = amgSystem.level(level).Ls\(f - amgSystem.level(level).U*u);
		u = amgSystem.level(level).Us\(f - amgSystem.level(level).L*u);
	end

end


end