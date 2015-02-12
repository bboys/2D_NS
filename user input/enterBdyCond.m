function [feMesh] = enterBdyCond(feMesh)
% adds information to feMesh.boundary describing for each part of the boundary
% its type and if it is of Dirichlet type, it contains a function giving the
% value along that part of the boundary
nrParts = size(feMesh.boundary,2);
for bdyPart = 1:nrParts
	feMesh.boundary(bdyPart).type = default({...
		'Choose the type of boundary condition', 'Dirichlet', 'Neumann'}, 1);
	if feMesh.boundary(bdyPart).type == 1
		% ask user to enter name of a predefined function (2 inputs 2 outputs)
		% or if default value is used => homogeneous dirichlet (no change)
		funcName = default('Name of Dirichlet function', 0, '', 'string');
		if strcmp(class(funcName), 'char') == 1
			feMesh.boundary(bdyPart).func =...
				str2func(['@(x,y)', funcName, '(x,y)']);
		else
			feMesh.boundary(bdyPart).func = [0; 0]; % homogeneous dirichlet
		end

	end
end

end