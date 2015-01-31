 function [] = plotMesh(feMesh)

figure
basisOrder = feMesh.basisOrder;
if strcmp(feMesh.meshType, 'quad')

	% only plot corners
	cornerIdx = [1, basisOrder + 1, (basisOrder + 1)^2,...
		(basisOrder + 1)*(basisOrder) + 1]
	patch('vertices', feMesh.node', 'faces', feMesh.elt(cornerIdx, :)',...
		'facecolor', 'w', 'linewidth', 1)


	% plot all nodes
	hold on
	plot(feMesh.node(1,:), feMesh.node(2,:), 'ks', 'MarkerSize', 2, ...
		'MarkerFaceColor', 'k')
	hold off
	
	% build title
	nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
	basisStr = ['Q',num2str(feMesh.basisOrder),...
		'P',num2str(feMesh.basisOrder - 1)];
	title(sprintf(['Quadrilateralisation by %d elements using ',basisStr,...
	 ' basis functions'], nrElts))
end
