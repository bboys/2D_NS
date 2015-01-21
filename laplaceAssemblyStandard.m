function [ D ] = laplaceAssemblyStandard(feMesh, localStiff)
% outputs the "laplace" matrix D
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
D = sparse(2*nrNodes,2*nrNodes);
for elt = 1:(feMesh.problemSize(1)*feMesh.problemSize(2))

    localNodes = feMesh.elt(:,elt);
    localDiff = localStiff.xx/feMesh.eltSize(1,elt)^2 +...
     localStiff.yy/feMesh.eltSize(2,elt)^2;

    D(localNodes,localNodes) = D(localNodes,localNodes) +...
     feMesh.area(elt)*localDiff;

    localNodes = localNodes + nrNodes;
    
    D(localNodes,localNodes) = D(localNodes,localNodes) +...
     feMesh.area(elt)*localDiff;

    % this simple representation works only with affine interpolation
    % equivalent quadrilateralations (definition of mesh).
end

end