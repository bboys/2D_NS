function [vorticity, vortMeshX, vortMeshY] = calcVorticity(feMesh, solVec)
% calculates the vorticity dv/dx - du/dy (not using basis functions)
% at the centers
nX = feMesh.problemSize(3);
nY = feMesh.problemSize(4);
nrNodes = nX*nY;
u = reshape(solVec(1:nrNodes), nY, nX);
v = reshape(solVec(nrNodes + 1:2*nrNodes), nY, nX);
dX = feMesh.gridX(2:end) - feMesh.gridX(1:end - 1); 
dY = feMesh.gridY(2:end) - feMesh.gridY(1:end - 1);
dUdY =  (u(2:end, 1:end - 1) - u(1:end - 1, 1:end - 1) + ...
	u(2:end, 2:end) - u(1:end - 1, 2:end))./(2*repmat(dY', 1, nX - 1));
dVdX =  (v(1:end - 1, 2:end) - v(1:end - 1, 1:end - 1) + ...
	v(2:end, 2:end) - v(2:end, 1:end - 1))./(2*repmat(dX, nY - 1, 1));
vorticity = dVdX - dUdY;

vortMeshX = (feMesh.gridX(1:end - 1) + feMesh.gridX(2:end))/2;
vortMeshY = (feMesh.gridY(1:end - 1) + feMesh.gridY(2:end))/2;
end