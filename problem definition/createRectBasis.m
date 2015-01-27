function [ localMatrix, basisOrder ] = createRectBasis(basisOrder )
%CREATERECTBASIS Summary of this function goes here
%   Detailed explanation goes here
% create matrices (if they don't exist) 

if nargin == 0
	basisQuery = {'Order of basis functions to use','Q1P0','Q2P1'};
	basisComments = {'stabilised', ''};
	basisOrder = default(basisQuery, 2, basisComments);
end

if basisOrder == 2
	if ~exist('localMatricesQ2P1.mat')
		localMatrix = createBasisQ2P1();
		save('localMatricesQ2P1.mat','localMatrix');
	else
		load localMatricesQ2P1.mat
	end
elseif basisOrder == 1
	if ~exist('localMatricesQ1P0.mat')
		localMatrix = createBasisQ1P0();
		save('localMatricesQ1P0.mat','localMatrix');
	else
		load localMatricesQ1P0.mat
	end
end


end

