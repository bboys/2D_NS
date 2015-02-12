function [] = cavityCenterPlot(solVec, feMesh, stokesSol, Re, lidVel)

% plot stuff

nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

solArray = reshape(solVec(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerSol = solArray(:,round(feMesh.problemSize(3)/2));

stokesArray = reshape(stokesSol(1:nrNodes),feMesh.problemSize(4),feMesh.problemSize(3));
centerStokes = stokesArray(:,round(feMesh.problemSize(3)/2));

% plotSol(feMesh, solVec, ['Re = ',num2str(Re)]);

yPoints = feMesh.node(2,1:feMesh.problemSize(4));

plot(yPoints,centerSol/lidVel,'r')
hold on
plot(yPoints,centerStokes/lidVel,'k')

legendStr = {'NS','Stokes', 'Stokes reference'};
% reference solution
load referenceSols
refSolStokes(:,2) = -refSolStokes(:,2)*2 + 1;

plot(1/2 + comsolCavityNSRe1(:,1)/2, comsolCavityNSRe1(:,2), 'kx')

if Re == 400
	plot(refSolYVal, refSol400, 'rx' )
	legendStr = [legendStr, {'NS reference, Re = 400'}];
elseif Re == 1000
	plot(1/2 + comsolNSRe1000(:,1)/2,  comsolNSRe1000(:,2), 'r+' )
	legendStr = [legendStr, {'NS reference, Re = 1000'}];
end




title(['Centerline velocity, Re = ',num2str(Re)])
legend(legendStr)

hold off
% grid on

end