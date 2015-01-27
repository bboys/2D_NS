function [ D ] = diffusionAssembly(feMesh, localStiff)
% outputs the "diffusion" matrix D, symmetric using Sij
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);

nrVBasisF = size(localStiff.xx,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;

localRangeI = repmat(1:nrVBasisF,nrVBasisF,1)'; 
localRangeI = localRangeI(:);
localRangeJ = repmat(1:nrVBasisF,nrVBasisF,1); 
localRangeJ = localRangeJ(:);

localStiffxx = reshape(localStiff.xx', nrVBasisFS,1);
localStiffxy = reshape(localStiff.xy', nrVBasisFS,1);
localStiffyx = reshape(localStiff.yx', nrVBasisFS,1);
localStiffyy = reshape(localStiff.yy', nrVBasisFS,1);

I1 = feMesh.elt(localRangeI,:); I1 = I1(:);
J1 = feMesh.elt(localRangeJ,:); J1 = J1(:);

I = [I1; nrNodes + I1; I1; nrNodes + I1];
J = [J1; J1; nrNodes + J1; nrNodes + J1];

K1 = kron(localStiffxx, 2*feMesh.area'./feMesh.eltSize(1,:).^2) +...
	kron(localStiffyy, feMesh.area'./feMesh.eltSize(2,:).^2);
K2 = repmat(localStiffyx,1,nrElts);
K3 = repmat(localStiffxy,1,nrElts);	
K4 = kron(localStiffyy, 2*feMesh.area'./feMesh.eltSize(2,:).^2) +...
	kron(localStiffxx, feMesh.area'./feMesh.eltSize(1,:).^2);

K = [K1(:); K2(:); K3(:); K4(:)];

D = sparse(I, J, K, 2*nrNodes, 2*nrNodes); % /2;
end