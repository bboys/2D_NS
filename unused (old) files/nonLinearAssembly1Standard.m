function [ N1 ] = nonLinearAssembly1Standard(feMesh, localNonLin, u)
% creates global nonlinear matrix based upon u (in the gradient term).
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrVBasisF = size(localNonLin.x,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;
N1 = sparse(2*nrNodes, 2*nrNodes);
for elt = 1:(feMesh.problemSize(1)*feMesh.problemSize(2))
        localNodes = [feMesh.elt(:,elt); nrNodes + feMesh.elt(:,elt)];
        u1 = u(localNodes(1:nrVBasisF)); % horizontal velocity component
        u2 = u(localNodes(nrVBasisF + 1:2*nrVBasisF)); % vertical ..  
        N11 = product3Didx2(localNonLin.x,u1, nrVBasisF)*feMesh.eltSize(2,elt);
        N12 = product3Didx2(localNonLin.x,u2, nrVBasisF)*feMesh.eltSize(2,elt);
        N13 = product3Didx2(localNonLin.y,u1, nrVBasisF)*feMesh.eltSize(1,elt);
        N14 = product3Didx2(localNonLin.y,u2, nrVBasisF)*feMesh.eltSize(1,elt);
        localN1 = [N11,N12';N13',N14];
        N1(localNodes,localNodes) = N1(localNodes,localNodes) +...
        	localN1'; 

 
        % this simple representation works only with affine interpolation
        % equivalent quadrilateralations (definition of mesh).
end
end

function [ mtx ] = product3Didx2(array3D, v, nrVBasisF)
% performs product and summation over 2nd index
% mtx(k,l) = sum_{r=1}^nrVBasisF array3D(k,r,l) v(r) 
mtx = zeros(nrVBasisF,nrVBasisF);
for r = 1:nrVBasisF
	mtx = mtx + reshape(v(r) * array3D(:,r,:),nrVBasisF,nrVBasisF);
end
end