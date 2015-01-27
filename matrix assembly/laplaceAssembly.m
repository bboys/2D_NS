function [ D ] = laplaceAssembly(feMesh, localStiff)
% outputs the laplace "diffusion" matrix D
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

nrVBasisF = size(localStiff.xx,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;

localRangeI = repmat(1:nrVBasisF,nrVBasisF,1)'; 
localRangeI = localRangeI(:);
localRangeJ = repmat(1:nrVBasisF,nrVBasisF,1); 
localRangeJ = localRangeJ(:);

localStiffxx = reshape(localStiff.xx',nrVBasisFS,1);
localStiffyy = reshape(localStiff.yy',nrVBasisFS,1);

I1 = feMesh.elt(localRangeI,:); I1 = I1(:);
J1 = feMesh.elt(localRangeJ,:); J1 = J1(:);

I = [I1;nrNodes + I1];
J = [J1;nrNodes + J1];

K1 = kron(localStiffxx, feMesh.area'./feMesh.eltSize(1,:).^2) +...
	kron(localStiffyy, feMesh.area'./feMesh.eltSize(2,:).^2);

K = [K1(:); K1(:)];

D = sparse(I, J, K, 2*nrNodes, 2*nrNodes);
end