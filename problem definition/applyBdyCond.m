function [nodeType, solVec] = applyBdyCond(feMesh, basisType)
% based on feMesh.boundary, this function outputs a struct nodeType 
% which defines for each boundary node whether or not it is free and moreover it
% initializes the solution vector 'solVec' to satisfy (if any) Dirichlet bdy.
% conditions
% feMesh.boundary(i) contains 
% (1) edge list (feMesh.boundary(i).edges), 
% (2) boundary type (feMesh.boundary(i).type),
% (3) (if any) Dirichlet value (function of the two spatial variables, or a 
% constant)
% boundary.func should be of size 2 x 1 (for u and v)

% number of parts of the boundary (disjoint boundary parts)
nrParts = size(feMesh.boundary,2);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

freeNodes = [];
freeVel = [];
freePressure = [];

fixedVel = [];
fixedPressure = [];


% maybe use
% types of boundary conditions (see page 43, Segal)
% 1 : dirichlet, u = g1,
% 2 : normal velocity u_n = g2, sigma_nt = g3 = 0, 
% 3 : tangential velocity u_t = g4, sigma_nn = g5 = 0, 
% 4 : sigma_nt = 0 = sigma_nn 

if (strcmp(basisType, 'Crouzeix-Raviart') == 1)
	nrPBasisF = 1/2*(feMesh.basisOrder)*(feMesh.basisOrder + 1);

	% initialize
	solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
	for bdyPart = 1:nrParts
		localNodes = unique(feMesh.boundary(bdyPart).edges(:));
		switch feMesh.boundary(bdyPart).type
		case 1
			% dirichlet
			if strcmp(class(feMesh.boundary(bdyPart).func),'function_handle' )
				tempSol = feMesh.boundary(bdyPart).func(...
					feMesh.node(1,localNodes), feMesh.node(2,localNodes));
				solVec(localNodes) = tempSol(1,:); % velocity in x dir
				solVec(localNodes + nrNodes) = tempSol(2,:); % and y dir
			elseif strcmp(class(feMesh.boundary(bdyPart).func),'double' )
				solVec(localNodes) = feMesh.boundary(bdyPart).func(1); 
				solVec(localNodes + nrNodes) = feMesh.boundary(bdyPart).func(2); 
			end
			fixedVel = [fixedVel; localNodes; localNodes + nrNodes];
		case 2

		case 3

		case 4


		end


	end
	% 'corners' were counted double
	fixedVel = unique(fixedVel);


	% fix one pressure value to zero
	fixedPressure = [1];
	freePressure = [2:nrPBasisF*nrElts]';

	freeVel = setdiff([1:2*nrNodes]', fixedVel);


	freeSol = [freeVel; 2*nrNodes + freePressure]; 


end

nodeType = struct('freeVel', freeVel, 'freePressure', freePressure,...
	'freeSol', freeSol, 'fixedVel', fixedVel, 'fixedPressure', fixedPressure);

end