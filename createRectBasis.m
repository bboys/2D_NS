function [ localMatrix ] = createRectBasis(basisOrder )
%CREATERECTBASIS Summary of this function goes here
%   Detailed explanation goes here
% create matrices (if they don't exist) 
if basisOrder == 'Q2P1'
	if ~exist('localMatricesQ2P1.mat')
		localMatrix = createBasisQ2P1();
		save('localMatricesQ2P1.mat','localMatrix');
	else
		load localMatricesQ2P1.mat
	end
elseif basisOrder == 'Q1P0'
	if ~exist('localMatricesQ1P0.mat')
		localMatrix = createBasisQ1P0();
		save('localMatricesQ1P0.mat','localMatrix');
	else
		load localMatricesQ1P0.mat
	end
end


end

