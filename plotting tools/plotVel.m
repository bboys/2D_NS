function [] = plotVel(feMesh, solVec, figureTitle, plotOption)
% plots the velocity field (magnitude) 

nX = feMesh.problemSize(3);
nY = feMesh.problemSize(4);
nrNodes = nX*nY; % velocity nodes
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
basisOrder = feMesh.basisOrder;

vMagnitude = sqrt(solVec(1:nrNodes).^2 + solVec(nrNodes + 1:2*nrNodes).^2); 

if strcmp(plotOption, 'flat')
	% color per face
	% simply average over nodes per element
	colorPerElt = mean(vMagnitude(feMesh.elt),1)';
	% scale
	colorPerElt = 100*(colorPerElt - min(colorPerElt))/(max(colorPerElt) - min(colorPerElt));

	% extract corners of each element
	cornerIdx = [1, basisOrder + 1, (basisOrder + 1)^2,...
		(basisOrder + 1)*(basisOrder) + 1];
	p = patch('vertices', feMesh.node', 'faces', feMesh.elt(cornerIdx, :)',...
		'linewidth', 1, 'FaceColor', 'flat', 'FaceVertexCData',...
		colorPerElt, 'CDataMapping', 'direct');
elseif strcmp(plotOption, 'interp')

	% color per vertex (faces are interpolated) 
	% extract corners of each element
	cornerIdx = [1, basisOrder + 1, (basisOrder + 1)^2,...
		(basisOrder + 1)*(basisOrder) + 1];
	p = patch('vertices', feMesh.node', 'faces', feMesh.elt(cornerIdx, :)',...
		'linewidth', 1, 'FaceColor', 'interp', 'FaceVertexCData',...
		vMagnitude, 'CDataMapping', 'scaled');
end

colormap jet;
axis off;
title(['Velocity magnitude, ',figureTitle])


end