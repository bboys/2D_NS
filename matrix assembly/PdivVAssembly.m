function [ L ] = PdivVAssembly( feMesh, localPdivV)
% outputs the pressure mass matrix for quads
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrElt = (feMesh.problemSize(1)*feMesh.problemSize(2)); % equal to nr of elts

nrVBasisF = size(localPdivV.x,2); % determine number of velocity basis functions
nrPBasisF = size(localPdivV.x,1); % and for the pressure
nrBasisFP = nrVBasisF*nrPBasisF; % product

I1 = repmat([1:nrElt*nrPBasisF],nrVBasisF,1);
I1 = I1(:);
I = [I1; I1]; % pressure functions

localRange = repmat([1:nrVBasisF]', nrPBasisF, 1);
J1 = feMesh.elt(localRange(:),:); J1 = J1(:);
J = [J1; nrNodes + J1];  % velocity functions (same in u,v direction)

localPdivVx = reshape(localPdivV.x', nrBasisFP, 1);
localPdivVy = reshape(localPdivV.y', nrBasisFP, 1);
K1 = bsxfun(@times, localPdivVx, feMesh.eltSize(2,:));
K2 = bsxfun(@times, localPdivVy, feMesh.eltSize(1,:));
K = [K1(:); K2(:)];

L = sparse(I, J, K, nrPBasisF*nrElt, 2*nrNodes);
end
