function [] = plotSol(feMesh, solVec, figureTitle)
% plots the velocity field (streamlines) and the pressure given a solution
% soLVec on the locations defined by plotSteps 
nX = feMesh.problemSize(3);
nY = feMesh.problemSize(4);
nrNodes = nX*nY; % velocity nodes
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

figure
% velocity field plot
% plotSteps =  ceil([nX/32 nY/32]);
% reshapeIdx = reshape(1:nrNodes, nY, nX);
% plotIdx = reshapeIdx(1:plotSteps(2):end,1:plotSteps(1):end);
% plotIdx = plotIdx(:);

% quiver(feMesh.node(1,plotIdx)', feMesh.node(2,plotIdx)', solVec(plotIdx),...
%   	solVec(nrNodes + plotIdx))


% contourf(reshape(sqrt(solVec(1:nrNodes).^2 + solVec(nrNodes + 1:2*nrNodes).^2), feMesh.problemSize(4),...
% 	feMesh.problemSize(3)));

vMagnitude = reshape(sqrt(solVec(1:nrNodes).^2 + solVec(nrNodes + 1:2*nrNodes).^2), feMesh.problemSize(4),...
 	feMesh.problemSize(3));
% imagesc(vMagnitude(end:-1:1, :)); 
uimagesc(feMesh.gridX, feMesh.gridY, vMagnitude(end:-1:1, :))
colormap jet;
% colorbar; 
axis off;


% streamline(reshape(feMesh.node(1,:), nY, nX),...
% 	reshape(feMesh.node(2,:), nY, nX),...
% 	reshape(solVec(1:nrNodes), nY, nX),...
% 	reshape(solVec((nrNodes + 1):2*nrNodes), nY, nX),-1:0.1:1,0);

title(['Velocity magnitude, ',figureTitle])


% axis([min(feMesh.node(1,:)), max(feMesh.node(1,:)), ...
% 	min(feMesh.node(2,:)), max(feMesh.node(2,:))])

% pressure contour plot
% figure
% pressureVal = solVec(2*nrNodes + 1:2*nrNodes + nrElts);
% pressureVal = pressureVal - min(pressureVal);
% contourf(reshape(pressureVal,... 
%  	feMesh.problemSize(2), feMesh.problemSize(1)));
% colorbar;
% title('Pressure contours')

% vorticity plot
% figure
% [vorticity, vortMeshX, vortMeshY] = calcVorticity(feMesh, solVec);
% uimagesc(vortMeshX, vortMeshY, abs(vorticity(end:-1:1, :)))
% colormap jet;

end