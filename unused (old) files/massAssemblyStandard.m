function [ M ] = massAssembly( feMesh, localMass)
% outputs the mass matrix based on feMesh for Q2 elements (biquadratic).
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
M = sparse(2*nrNodes,2*nrNodes);
for elt = 1:(feMesh.problemSize(1)*feMesh.problemSize(2))
        localNodes = feMesh.elt(:,elt);
        M(localNodes,localNodes) = M(localNodes,localNodes) + ...
            feMesh.area(elt)*localMass;
        % this simple representation works only with affine interpolation
        % equivalent quadrilateralations (definition of mesh).
end
% same basis functions for u and v
M(nrNodes+1:end,nrNodes+1:end) = M(1:nrNodes,1:nrNodes);
end
