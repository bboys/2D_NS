function [ M ] = vmassAssembly( feMesh, localMass)
% outputs the mass matrix based on feMesh for Q2 elements (biquadratic).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

nrVBasisF = size(localMass,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;

localRange1 = repmat(1:nrVBasisF,nrVBasisF,1)'; 
localRange1 = localRange1(:);
localRange2 = repmat(1:nrVBasisF,nrVBasisF,1); 
localRange2 = localRange2(:);

I = feMesh.elt(localRange1,:);
J = feMesh.elt(localRange2,:);
localMassVec = reshape(localMass',nrVBasisFS,1);
K = kron(localMassVec, feMesh.area');
M1 = sparse(I(:),J(:),K(:),nrNodes,nrNodes); % M1(I(k),J(k)) = K(k)

% same basis functions for u and v
M = sparse(2*nrNodes,2*nrNodes);
M(1:nrNodes, 1:nrNodes) = M1;
M(nrNodes+1:end, nrNodes+1:end) = M1;

end
