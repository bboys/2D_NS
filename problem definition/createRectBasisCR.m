function [ localMatrix ] = createRectBasisCR(basisOrder)
% creates the local matrices for arbitrary basis order k Q(k)Pk - 1
% where k = basisOrder. 
% using Crouzeix-Raviart quadrilaterals
nrVBasisF = (basisOrder + 1)^2; % dim(Q_k)
nrPBasisF = 1/2*(basisOrder)*(basisOrder + 1); % dim(P_k-1)


% define reference element
points = linspace(0, 1, basisOrder + 1); % different choise?

phi = zeros(basisOrder + 1, basisOrder + 1); % each row contains coeff basisF
nodeFunc = [ones(1, basisOrder + 1); -points]; 
for i = 1:basisOrder + 1
	phi(i, :) = multiPolyProduct(nodeFunc(:,[1:i-1,i+1:basisOrder + 1]));
	phi(i, :) = phi(i, :)/polyval(phi(i, :), points(i));

	% % display  basis
	% x = linspace(0,1,200);
	% plot(x,polyval(phi(i, :), x))
	% hold on
	% pause
end

% local velocity mass matrix
% one dimensional mass matrix (since basisF(ij) = phii(x)*phij(y)) we can
% split integration in x and y direction
vMass1D = zeros(basisOrder + 1, basisOrder + 1);
for i = 1:basisOrder + 1
	for j = i:basisOrder + 1 % use symmetry
		basisF1D = polyProduct(phi(i,:), phi(j, :));
		vMass1D(i, j) = sum(polyint(basisF1D')); % integrate from 0 to 1...
	end
end
% duplicate due to symmetry
tempDiag = diag(vMass1D); 
vMass1D = vMass1D + vMass1D'; 
vMass1D(1:(basisOrder + 2):end) = tempDiag; 
localMatrix.vmass = kron(vMass1D, vMass1D);

% local velocity stiffness matrix
% derivatives of 1D basis functions
phix = zeros(basisOrder + 1, basisOrder);
for i = 1:basisOrder + 1
	phix(i, :) = polyder(phi(i, :));
end

vStiffxx1D = zeros(basisOrder + 1, basisOrder + 1); % int (d lambdai dx)*(d lambdaj dx)
vStiffx1D = zeros(basisOrder + 1, basisOrder + 1); % int lambdai (d lambdaj dx)
for i = 1:basisOrder + 1
	for j = 1:basisOrder + 1
		vStiffxx1D(i, j) = sum(polyint(polyProduct(phix(i,:), phix(j, :))'));
		vStiffx1D(i, j) = sum(polyint(polyProduct(phi(i,:), phix(j, :))'));
	end
end
localMatrix.stiff.xx = kron(vStiffxx1D, vMass1D);
localMatrix.stiff.yy = kron(vMass1D, vStiffxx1D);
localMatrix.stiff.xy = kron(vStiffx1D', vStiffx1D);
localMatrix.stiff.yx = kron(vStiffx1D, vStiffx1D');

% local nonlinear term
% lnonlin.x(k,r,l) = (int_0^1)^2 basisF(k) d/dx basisF(r) basisF(l) dx
% utilize tensor product: basisF(ij) = phii(x)*phij(y)

vNonLin1Dx = zeros(basisOrder + 1, basisOrder + 1, basisOrder + 1);
% vNonLin1Dx(i,j,k) = int phii * (d phij dx) * phik
vNonLin1D = zeros(basisOrder + 1, basisOrder + 1, basisOrder + 1);
% vNonLin1D(i,j,k) = int phii * phij * phik
for i = 1: basisOrder + 1
	for j = 1: basisOrder + 1
		for k = 1: basisOrder + 1
			tempProd = polyProduct(phi(i,:), phi(k,:));
			vNonLin1Dx(i,j,k) = sum(polyint(polyProduct(tempProd, phix(j,:))'));
			vNonLin1D(i,j,k) = sum(polyint(polyProduct(tempProd, phi(j,:))'));
		end
	end
end

% localNonLin.x(k,r,l) = (int_0^1)^2 basisF(k) d/dx basisF(r) basisF(l) dx
localMatrix.nonlin.x = superkron(vNonLin1Dx, vNonLin1D);
localMatrix.nonlin.y = superkron(vNonLin1D, vNonLin1Dx);

% pressure basis functions live in P_(basisOrder - 1), one nodal point, and 
% up to all the basisOrder (mixed) derivatives, all evaluated in center
% hence basis functions are functions as they appear in 2D taylor series (around 1/2 1/2)
% that is, Pij = (x-1/2)^i*(y-1/2)^j, integrals can again be split => 1D
% for 0 <= i + j <= basisOrder

% let lambda_l =  (x-1/2)^l, since then
% int Pij Pmn dxdy = (int lambdai lambdam dx)*(int lambdaj lambdan dx)/((i+j)!(m+n)!)

% hence we compute int_01 (x-1/2)^l, for l = 0, ... ,2*basisOrder
% note that int_01 (x-1/2)^l = (2^(-l-1))/(l+1) * (1 + (-1)^(-l))
l = 0:2*(basisOrder - 1);
pMass1Dvec = ((2.^(-l-1))./(l+1)).*(1 + (-1).^(-l));

% now construct local matrix pMass1D_ij = int_01 lambdai lambdaj dx, that is
% pMass1D_(i+1)(j+1) = pMass1Dvec(i+j), i,j=0,...,basisOrder
pMass = zeros(nrPBasisF, nrPBasisF);
I = 0;
for i = 0:basisOrder-1
	for j = 0:basisOrder-1-i
		I = I + 1;
		J = 0;
		for m = 0:basisOrder-1
			for n = 0:basisOrder-1-m
				J = J + 1;
				pMass(I,J) = pMass1Dvec(i+m+1)*pMass1Dvec(j+n+1)/...
					(factorial(i+j)*factorial(m+n));

			end
		end
	end
end % vectorize?
localMatrix.pmass = pMass;

% finally localPdivV.x_r(i,j) = int_0^1 int_0^1  d/dx_r basisF(i) * basisP(j)
% we will write the pressure 1D basis functions in terms of the velocity basis
% functions (standard basis 1, x, x^2 ...). hence we calculate the change of basis
temp = 0:basisOrder - 1;
changeBasis = abs(pascal(basisOrder ,1)).*((1/2).^bsxfun(@minus,temp', temp ));

changeBasis = inv(changeBasis);
% now each row i of changeBasis is representation of lambdai in terms of 
% standard basis

PdivV1D = zeros(basisOrder, basisOrder + 1);
PdivV1Dx = zeros(basisOrder, basisOrder + 1);

for i = 1:basisOrder
	for j = 1:basisOrder + 1
		pFunc = changeBasis(i, end:-1:1);

		PdivV1D(i, j) = diff(polyval(polyint(polyProduct(phi(j, :), pFunc)'),[0 1]));
		PdivV1Dx(i, j) = diff(polyval(polyint(polyProduct(phix(j, :), pFunc)'), [0,1]));
	end
end


localMatrix.pdivv.x = zeros(nrPBasisF,nrVBasisF);
localMatrix.pdivv.y = zeros(nrPBasisF,nrVBasisF);

I = 0;
for i = 0:basisOrder-1
	for j = 0:basisOrder-1-i
		I = I + 1;
		for k = 1:basisOrder + 1 % loop over velocity basis functions
			for l = 1:basisOrder + 1
				localMatrix.pdivv.x(I, (k-1)*(basisOrder + 1) + l) = ...
					PdivV1D(i+1, l)*PdivV1Dx(j+1, k)/factorial(i+j);

				localMatrix.pdivv.y(I, (k-1)*(basisOrder + 1) + l) = ...
					PdivV1Dx(i+1, l)*PdivV1D(j+1, k)/factorial(i+j);
			end
		end
	end
end % vectorize?

localMatrix.basisType = 'Crouzeix-Raviart';

end

function [P] = polyProduct(p1, p2)
% performs the product of two polynomials, given by p1, p2
deg1 = length(p1) - 1; deg2 = length(p2) - 1; 
p1 = p1(:); p2 = p2(:); % columnn vectors

P = zeros(deg1 + deg2 + 1, 1); % resulting product
for coef = 1:deg1 + 1
	P = P + p1(coef)*[zeros(coef - 1, 1); p2; zeros(deg1 - coef + 1, 1)];
end


end

function [P] = multiPolyProduct(Parray)
% performs the product of any collection of polynomials
P = Parray(:, 1);
for i = 2:size(Parray, 2)
	P = polyProduct(P, Parray(:, i));
end

end

