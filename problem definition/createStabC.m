function [stabC] = createStabC(feMesh, localMatrix)
stabC = [];

% for now only for the CR family of rect elements

if (strcmp(localMatrix.basisType,'Crouzeix-Raviart') == 1) 
	nX = feMesh.problemSize(1);
	nY = feMesh.problemSize(2);
	if (feMesh.basisOrder == 1)
		% subdivision is by 2x2 patches
		xDif = feMesh.gridX(3:2:end) - feMesh.gridX(1:2:end-2);
		yDif = feMesh.gridY(3:2:end) - feMesh.gridY(1:2:end-2);

		patchSize = bsxfun(@times, yDif, xDif');

		% create stabilisation matrix stabC
		beta = 1/4;
		stabPatch = [2 -1 0 -1;
					-1 2 -1 0; 
					0 -1 2 -1; 
					-1 0 -1 2];

		stabC = sparse(1:nX*nY/4, 1:nX*nY/4, patchSize(:)/4);
		stabC = beta*kron(stabC, stabPatch);
	else
		nrPBasisF = 1/2*(feMesh.basisOrder)*(feMesh.basisOrder + 1);
		stabC = sparse(nrPBasisF*nX*nY, nrPBasisF*nX*nY);
	end

end

end
