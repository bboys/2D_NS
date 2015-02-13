function [amgSystem] = createAmgSystem(A, levels, setup)
% creates struct amgSystem containing the matrices 
% A^h, ..., A^(2^levels h), and corresponding interpolation
% operators
%
% amgSystem.level(i).matrix, and
% amgSystem.level(i).interp, such that
%
% I^{2h}_h = amgSystem.level(1).interp'
% amgSystem.level(1).matrix = A^h = A (original matrix, finest grid)
%
% levels indicates the maximum level (1 is sufficient for two cycle AMG)

% initialize
amgSystem.level(1).matrix = A;

% apply same setup on each level
for level = 1:levels
	amgSystem.level(level).interp = createInterpOp(amgSystem.level(level).matrix, setup);
	amgSystem.level(level + 1).matrix =...
		amgSystem.level(level).interp'*amgSystem.level(level).matrix*...
		amgSystem.level(level).interp;
end

end