function [ N1 ] = nonLinearAssembly1(feMesh, localNonLin, u)
% creates global nonlinear matrix based upon u (in the gradient term).
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

I = [I1; I1; I1 + nrNodes; I1 + nrNodes];
J = [J1; J1 + nrNodes; J1; J1 + nrNodes];

K1 = zeros(nrVBasisFS*nrElts,1);
K2 = zeros(nrVBasisFS*nrElts,1);
K3 = zeros(nrVBasisFS*nrElts,1);
K4 = zeros(nrVBasisFS*nrElts,1);


for r = 1:nrVBasisF
	localNonLinx = reshape(localNonLin.x(:,r,:),nrVBasisF,nrVBasisF);
	localNonLiny = reshape(localNonLin.y(:,r,:),nrVBasisF,nrVBasisF);

	% localNonLinx = localNonLinx(:);
	% localNonLiny = localNonLiny(:);

	% K1temp = bsxfun(@times, localNonLinx, [feMesh.eltSize(2,:).*...
	% 	u1(feMesh.elt(r,:))', feMesh.eltSize(2,:).*u2(feMesh.elt(r,:))']); 
	% K1temp = K1temp(:);



	% K4temp = bsxfun(@times, localNonLiny, [feMesh.eltSize(1,:).*...
	% 	u1(feMesh.elt(r,:))', feMesh.eltSize(1,:).*u2(feMesh.elt(r,:))']); 
	% K4temp = K4temp(:);

	K1temp = bsxfun(@times, localNonLinx(:), feMesh.eltSize(2,:).*...
		u1(feMesh.elt(r,:))'); 
	K1temp = K1temp(:);

	localNonLinx = localNonLinx';

	K2temp = bsxfun(@times, localNonLinx(:), feMesh.eltSize(2,:).*...
		u2(feMesh.elt(r,:))'); 
	K2temp = K2temp(:);


	K3temp = bsxfun(@times, localNonLiny(:), feMesh.eltSize(1,:).*...
		u1(feMesh.elt(r,:))'); 
	K3temp = K3temp(:);

	localNonLiny = localNonLiny';

	K4temp = bsxfun(@times, localNonLiny(:), feMesh.eltSize(1,:).*...
		u2(feMesh.elt(r,:))'); 
	K4temp = K4temp(:);


	K1 = K1 + K1temp;
	K2 = K2 + K2temp;
	K3 = K3 + K3temp;
	K4 = K4 + K4temp;

end

K = [K1; K3; K2; K4];

% combine four blocks
N1 = sparse(I,J,K,2*nrNodes,2*nrNodes);

end