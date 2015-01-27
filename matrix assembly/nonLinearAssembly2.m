function [ N2 ] = nonLinearAssembly2(feMesh, localNonLin, u)
% creates global nonlinear matrix based upon u (in the first term).
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4); % velocity nodes
nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
u1 = u(1:nrNodes); u2 = u((1 + nrNodes):2*nrNodes);

nrVBasisF = size(localNonLin.x,1); % determine number of velocity basis functions
nrVBasisFS = nrVBasisF^2;

localRangeI = repmat(1:nrVBasisF,nrVBasisF,1)'; % 123123123
localRangeI = localRangeI(:);
localRangeJ = repmat(1:nrVBasisF,nrVBasisF,1); % 111222333
localRangeJ = localRangeJ(:);

I1 = feMesh.elt(localRangeI,:); I1 = I1(:);
J1 = feMesh.elt(localRangeJ,:); J1 = J1(:);

I = [I1; I1 + nrNodes];
J = [J1; J1 + nrNodes];

K = zeros(nrVBasisFS*nrElts,1);

% bsxfun could be usefull here !!!!

for r = 1:nrVBasisF
        localNonLinx = reshape(localNonLin.x(:,:,r),nrVBasisFS,1);
        localNonLiny = reshape(localNonLin.y(:,:,r),nrVBasisFS,1);

        K1temp = bsxfun(@times, localNonLinx, feMesh.eltSize(2,:).*...
                u1(feMesh.elt(r,:))'); 
        K1temp = K1temp(:);

        K2temp = bsxfun(@times, localNonLiny, feMesh.eltSize(1,:).*...
                u2(feMesh.elt(r,:))'); 
        K2temp = K2temp(:);

        K = K + K1temp + K2temp;
end

K = [K; K]; % should be divided by two, but why? (or diffusion matrix *2)
N2 = sparse(I,J,K,2*nrNodes,2*nrNodes);
end