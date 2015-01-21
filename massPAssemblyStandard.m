function [ L ] = massPAssemblyStandard( feMesh, localPMass)
% outputs the " pressire mass" matrix relating the pressure (P1) and velocity 
% (Q2) unknowns
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrPNodes = (feMesh.problemSize(1)*feMesh.problemSize(2)); % equal to nr of elts
L = sparse(3*nrPNodes, 2*nrNodes);
for elt = 1:nrPNodes
        localNodes = feMesh.elt(:,elt);
        pNodes = elt + [0,nrPNodes,2*nrPNodes]; % pressure node (3x)
        L(pNodes,localNodes) = ...
         feMesh.area(elt)*localPMass.x/feMesh.eltSize(1,elt);
        L(pNodes,localNodes + nrNodes) = ...
         feMesh.area(elt)*localPMass.y/feMesh.eltSize(2,elt);
        % this simple representation works only with affine interpolation
        % equivalent quadrilateralations (definition of mesh).
end

end