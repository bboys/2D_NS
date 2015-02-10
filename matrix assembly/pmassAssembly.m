function [ Q ] = pmassAssembly( feMesh, localMass)
% outputs the mass matrix based on feMesh for Q2 elements (biquadratic).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);
nrElt = feMesh.problemSize(1)*feMesh.problemSize(2);

nrPBasisF = size(localMass,1);
nrPBasisFS = nrPBasisF^2;


I = repmat([1:nrElt*nrPBasisF],nrPBasisF,1);
I = reshape(I, nrElt*nrPBasisF,nrPBasisF);
I = I(:);

J = repmat(1:nrPBasisF, 1, nrPBasisF);
J = bsxfun(@plus, [0:nrPBasisF:(nrElt - 1)*nrPBasisF]', J)';
J = J(:);

localMassVec = reshape(localMass',nrPBasisFS,1);
K = bsxfun(@times, localMassVec, feMesh.area');

Q = sparse(I(:),J(:),K(:),nrElt*nrPBasisF,nrElt*nrPBasisF); 

% this function treats order of pressure basis functions per element, such that
% a nice structure is obtained
end
