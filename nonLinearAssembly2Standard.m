function [ N2 ] = nonLinearAssembly2Standard(feMesh, localNonLin, u)
% creates global nonlinear matrix based upon u (in the first term).
% basic implementation = slow (needs improvement).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrVBasisF = size(localNonLin.x,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;

N2 = sparse(2*nrNodes, 2*nrNodes);
for elt = 1:(feMesh.problemSize(1)*feMesh.problemSize(2))
        localNodes = [feMesh.elt(:,elt); nrNodes + feMesh.elt(:,elt)];
        u1 = u(localNodes(1:nrVBasisF)); % horizontal velocity component
        u2 = u(localNodes(nrVBasisF + 1:2*nrVBasisF)); % vertical ..  

        N21 = product3Didx1(localNonLin.x,u1, nrVBasisF)*feMesh.eltSize(2,elt) +...
                product3Didx1(localNonLin.y,u2, nrVBasisF)*feMesh.eltSize(1,elt);
        localN2 = [N21, zeros(nrVBasisF,nrVBasisF);...
                zeros(nrVBasisF,nrVBasisF),N21];

        N2(localNodes,localNodes) = N2(localNodes,localNodes) +...
                localN2'; 
        % this simple representation works only with affine interpolation
        % equivalent quadrilateralations (definition of mesh).
end
end

function [ mtx ] = product3Didx1(array3D, v, nrVBasisF)
% performs product and summation over 1st index
% mtx(r,l) = sum_{k=1}^nrVBasisF array3D(k,r,l) v(k) 
mtx = zeros(nrVBasisF,nrVBasisF);
for k = 1:nrVBasisF
        mtx = mtx + reshape(v(k) * array3D(k,:,:),nrVBasisF,nrVBasisF);
end

end