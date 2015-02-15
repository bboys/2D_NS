function [interpOp] = createInterpOp(Ah, setup)
% creates interpolation operator I^h_{2h} mapping from the coarse
% grid corresponding to A^{2h} to the fine grid
% setup.amg.coarseMethod defines which coarsening method is to be used
% setup.amg.interpMethod ... interpolation ...
% setup.amg.theta is the connectivity threshold

% temporary setup for testing
% Ah = gallery('poisson', 5);
% setup.amg.coarseMethod = 'RS';
% setup.amg.interpMethod = 'classical';
% setup.amg.theta = 0.6;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rng('default') 	  % temporary % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxiliary strength matrix, u_i strongly depends on u_j (u_j strongly influences u_i)
nrVar = size(Ah,1);
noDiagAh = Ah + spdiags(Inf*ones(nrVar, 1), 0, nrVar, nrVar);

rowMax = max(-noDiagAh, [], 2);
auxStrength = bsxfun(@ge, -Ah, setup.amg.theta*rowMax);
measure = full(sum(auxStrength)); % for selecting initial fine nodes
measure = measure(:);

% initial c, f and rest nodes
fNodes = measure == 0;
cNodes = logical(zeros(nrVar,1));
restNodes = logical(ones(nrVar,1));
restNodes(fNodes) = 0;

tempAuxStrength = auxStrength;
if setup.amg.coarseMethod == 1 
	% classical Ruge-Stuben coarsening (satisfies F-F condition)

	% first pass
	while any(restNodes)
		% select node attaining maximum
		[~, cRestNodes] = max(measure);

		% select 'strongly dependant' neighbours of cRestNodes (from restNodes)
		[~, fRestNodes] = find(tempAuxStrength(cRestNodes, :));

		% increase measure of nodes that strongly influence fRestNodes
		[~, strongFRest] = find(tempAuxStrength(fRestNodes, :));

		tempRemove = [cRestNodes(:); fRestNodes(:)];
		measure(tempRemove) = 0;
		tempAuxStrength(tempRemove, tempRemove) = 0;
		measure(strongFRest) = measure(strongFRest) + 1;

		% update node sets
		cNodes(cRestNodes) = 1;
		fNodes(fRestNodes) = 1;
		restNodes(tempRemove) = 0;
	end

	% second pass (enforce F-F condition)
	% connFF = auxStrength(fNodes, fNodes);

	% for each off-diagonal nonzero element, we must check if the corresponding
	% F nodes share a common C node.


elseif setup.amg.coarseMethod == 2
	% standard PMIS coarsening (does not guarrantee F-F condition)

	measure = measure + rand(nrVar, 1);

	% unite both types of neighbours (step may be skipped if symmetric Ah)
	tempAuxStrength = auxStrength | auxStrength';

	% give each node its measure
	strength =  bsxfun(@times, measure', tempAuxStrength);
	tempLTN = [1:nrVar]'; % convert logical to normal index
	while any(restNodes)

		% select those nodes which have measure larger than all neighbours
		cRestNodes = all(bsxfun(@gt, measure(restNodes),...
			strength(restNodes, restNodes)), 2);
		
		% convert indices
		restNodesN = tempLTN(restNodes);
		cRestNodes = restNodesN(cRestNodes);

		% select 'strongly dependant' neighbours of cRestNodes (from restNodes)
		[~, fRestNodes] = find(strength(cRestNodes, :));

		% update arrays
		tempRemove = [cRestNodes(:); fRestNodes(:)];

		measure(tempRemove) = 0;
		strength(tempRemove, tempRemove) = 0;

		% update node sets
		cNodes(cRestNodes) = 1;
		fNodes(fRestNodes) = 1;
		restNodes(tempRemove) = 0;
	end



elseif setup.amg.coarseMethod == 3
	% CLJP coarsening
end

% count
nrcNodes = sum(cNodes);

% the coarse interpolatory set: ith row gives the neighbours which strongly
% influence i AND are in the cNodes set (C_i)
coarseInterp = bsxfun(@and, auxStrength, cNodes');

% fNodeStrong: ith row gives neighbours which strongly influence i but are not 
% in the cNodes set (hence they are in the fNodes set) (D^s_i)
fNodeStrong = bsxfun(@and, auxStrength, fNodes');

% weakNeighbours: ith row gives neighbours which do not strongly influence i
% (the left over neighbours) (D^w_i)

zeroDiagAh = Ah - spdiags(diag(Ah), 0, nrVar, nrVar);
weakNeighbours = xor(logical(zeroDiagAh), (coarseInterp | fNodeStrong));

% general interpolation rule:
% (I^h_{2h}e)_i = e_i 							if i in cNodes, 
% 				= sum_{j in C_i} w_{ij} e_j 	if i in fNodes
% here C_i is defined by coarseInterp(i, :)

if setup.amg.interpMethod == 2
	% F-F interpolation (extend C_i and D^s_i)

	%  the strong F-F connections
	strongFF = sparse(nrVar, nrVar);
	strongFF(fNodes, fNodes) = fNodeStrong(fNodes, fNodes);

	% if strongFF(i,j) is nonzero, then we have an fNode strongly influencing
	% an fNode, we must check whether they share a cNode which strongly influences
	% both i and j. If such a cNode does not exist, then we must add C_j to C_i
	% This occurs if (coarseInterp(:,i) & coarseInterp(:,j)) is all zeros
	[ffRows, ffCols] = find(strongFF);

	coarseI = coarseInterp(ffRows, :); % per ffNode this gives all strongly influencing cNodes
	coarseJ = coarseInterp(ffCols, :);

	% if badFF(l) == 1, then ffRows(l) and ffCols(l) have a strong F-F connection
	% without sharing a cNode
	badFF = full(~any(coarseI & coarseJ, 2));
	
	% extend coarseInterp (add C_j to C_i)
	coarseInterp(ffRows(badFF), :) = coarseInterp(ffRows(badFF), :) |...
		coarseInterp(ffCols(badFF), :);

end
% if not F-F then classical (skip extension of sets)


if setup.amg.interpMethod == 1 | setup.amg.interpMethod == 2


	% requires some attention (too much is calculated now, and divided by zero)
	% storage of NaN?
	
	% precompute denominator: denom(i) = aii + sum_{n in D^w_i} a_{in}, 
	% for i = 1:nrVar, (result is a column vector)
	denom = diag(Ah) + sum(Ah.*weakNeighbours,2);


	% precompute denominator in numerator: denomNum(i,m) = sum_{k in C_i} a_{mk},
	% for i,m = 1:nrVar (result is a matrix)
	denomNum = coarseInterp*Ah';

	% compute numerator: 
	% numer(i,j) = a_{ij} + sum_{m in D^s_i (a_im a_mj)/denomNum(i,m)} 

	numer = Ah;
	temp1 = sparse(nrVar, nrVar);
	temp1(fNodeStrong) = Ah(fNodeStrong)./denomNum(fNodeStrong);

	temp1(isnan(temp1)) = 0;
	temp1(isinf(temp1)) = 0;

	numer = numer + temp1*Ah;

	% numer = Ah;
	% for i = 1:nrVar
	% 	for m = 1:nrVar
	% 		if fNodeStrong(i, m) == 1 & Ah(i, m) ~= 0
	% 			for j = 1:nrVar
	% 				if Ah(m, j) ~= 0
	% 					numer(i,j) = numer(i,j) + Ah(i, m)*Ah(m, j)/denomNum(i,m);
	% 				end
	% 			end
	% 		end
	% 	end
	% end


	% resulting weights: interpOp(i,j) = - (numer(i,j)/denom(i))
	interpOp = - bsxfun(@times, numer(:, cNodes), 1./denom);

	% coarse node is interpolated by itself
	interpOp(cNodes, :) =  speye(nrcNodes);

end



end