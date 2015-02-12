function [ localMatrix, basisOrder, basisType ] = createBasis(basisOrder)
% ask user to give the order of the basis functions to be used

if nargin == 0
	% basisQuery = {'Order of basis functions to use','Q1P0','Q2P1','Q3P2','Q4P3', 'Q5P4'};
	% basisComments = {'stabilised', '', '', '', ''};
	% basisOrder = default(basisQuery, 2, basisComments);
	basisOrder = default('Order of basis functions to use Q(k)P(k-1)', 2);

end

basisType = 'Crouzeix-Raviart';
[ localMatrix ] = createRectBasisCR(basisOrder);

end
