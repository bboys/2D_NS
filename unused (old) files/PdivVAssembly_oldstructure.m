function [ L ] = PdivVAssembly_oldstructure( feMesh, localPdivV)
% outputs the " pressure mass" matrix relating the pressure (P1) and velocity 
% (Q2) unknowns
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrPNodes = (feMesh.problemSize(1)*feMesh.problemSize(2)); % equal to nr of elts

nrVBasisF = size(localPdivV.x,2); % determine number of velocity basis functions
nrPBasisF = size(localPdivV.x,1); % and for the pressure
nrBasisFP = nrVBasisF*nrPBasisF;


I1 = reshape([1:nrPBasisF*nrPNodes],nrPNodes,nrPBasisF)';
I1 = repmat(I1(:),1,nrVBasisF)';
I1 = I1(:);

I = [I1; I1];


localRange2 = repmat(1:nrVBasisF,1,nrPBasisF);

J1 = feMesh.elt(localRange2,:); J1 = J1(:);
J = [J1; nrNodes + J1];

localPdivVx = reshape(localPdivV.x', nrBasisFP, 1);
localPdivVy = reshape(localPdivV.y', nrBasisFP, 1);
K1 = kron(localPdivVx, feMesh.eltSize(2,:));
K2 = kron(localPdivVy, feMesh.eltSize(1,:));
K = [K1(:); K2(:)];

L = sparse(I, J, K, nrPBasisF*nrPNodes, 2*nrNodes);
end