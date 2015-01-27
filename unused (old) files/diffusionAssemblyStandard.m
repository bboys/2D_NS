function [ D ] = diffusionAssemblyStandard(feMesh, localStiff)
% outputs the "diffusion" matrix D for Q2 elements
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
D = sparse(2*nrNodes,2*nrNodes);
for elt = 1:(feMesh.problemSize(1)*feMesh.problemSize(2))
	% include cross terms
    localNodes = [feMesh.elt(:,elt),feMesh.elt(:,elt) + nrNodes];
    L1 = 2*localStiff.xx/feMesh.eltSize(1,elt)^2 +...
     localStiff.yy/(feMesh.eltSize(2,elt)^2);
    L2 = localStiff.yx/(feMesh.eltSize(1,elt)*feMesh.eltSize(2,elt));
    L3 = localStiff.xy/(feMesh.eltSize(1,elt)*feMesh.eltSize(2,elt));
	L4 = 2*localStiff.yy/feMesh.eltSize(2,elt)^2 +...
     localStiff.xx/(feMesh.eltSize(1,elt)^2);
    localDiff = [L1,L2;L3,L4];
    D(localNodes,localNodes) = D(localNodes,localNodes) +...
     feMesh.area(elt)*localDiff;
    % this simple representation works only with affine interpolation
    % equivalent quadrilateralations (definition of mesh).
end

end