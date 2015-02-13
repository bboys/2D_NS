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

% auxiliary strength matrix, u_i strongly depends on u_j (u_j strongly influences u_i)
nrVar = size(Ah,1);
noDiagAh = Ah + spdiags(Inf*ones(nrVar, 1), 0, nrVar, nrVar);

rowMax = max(-noDiagAh, [], 2);
auxStrength = bsxfun(@ge, -Ah, setup.amg.theta*rowMax);
measure = full(sum(auxStrength)); % for selecting initial fine nodes

% initial c, f and rest nodes
fNodes = find(measure == 0);
fNodes = fNodes(:);
cNodes = [];
restNodes = setdiff(1:nrVar, fNodes);
if strcmp(setup.amg.coarseMethod, 'RS')
	% classical Ruge-Stuben coarsening (satisfies F-F condition)

	% first pass
	while ~isempty(restNodes)
		restStrength = auxStrength(restNodes, restNodes);
		restMeasure = measure(restNodes)';

		% select node attaining maximum
		[~, cRestNodes] = max(restMeasure);

		% select 'strongly dependant' neighbours of cRestNodes (from restNodes)
		[~, fRestNodes] = find(restStrength(cRestNodes, :));

		% increase measure of nodes that strongly influence fRestNodes
		[~, strongFRest] = find(restStrength(fRestNodes, :));

		% nodes in terms of position in restnodes
		strongFRest = restNodes(strongFRest);
		cRestNodes = restNodes(cRestNodes);
		fRestNodes = restNodes(fRestNodes);

		measure(strongFRest) = measure(strongFRest) + 1;

		% update node sets
		cNodes = [cNodes; cRestNodes(:)];
		fNodes = [fNodes; fRestNodes(:)];
		restNodes = setdiff(restNodes, [cRestNodes(:); fRestNodes(:)]);
	end

	% second pass (enforce F-F condition)
	connFF = auxStrength(fNodes, fNodes);

	% for each off-diagonal nonzero element, we must check if the corresponding
	% F nodes share a common C node.


elseif strcmp(setup.amg.coarseMethod, 'PMIS')
	% standard PMIS coarsening (does not guarrantee F-F condition)

	measure = measure + rand(1, nrVar);

	while ~isempty(restNodes)
		% select subset of restNodes of nodes which have a measure stronger than
		% both its strongly dependend neighbours and neighbours which strongly
		% depend on it
		restStrength = auxStrength(restNodes, restNodes);

		% unite both types of neighbours (step may be skipped if symmetric Ah)
		restStrength = restStrength | restStrength';

		% measure of restNodes (column array)
		restMeasure = measure(restNodes)';

		% select those nodes which have measure larger than all neighbours
		restStrength = bsxfun(@times, restMeasure', restStrength);
		cRestNodes = all(bsxfun(@gt, restMeasure, restStrength), 2);

		% select 'strongly dependant' neighbours of cRestNodes (from restNodes)
		[~, fRestNodes] = find(restStrength(cRestNodes, :));

		% nodes in terms of position in restnodes
		cRestNodes = restNodes(cRestNodes);
		fRestNodes = restNodes(fRestNodes);

		% update node sets
		cNodes = [cNodes; cRestNodes(:)];
		fNodes = [fNodes; fRestNodes(:)];
		restNodes = setdiff(restNodes, [cRestNodes(:); fRestNodes(:)]);
	end
elseif strcmp(setup.amg.coarseMethod, 'CLJP')
	% CLJP coarsening
end

% sort
cNodes = sort(cNodes);
fNodes = sort(fNodes);

% the coarse interpolatory set: ith row gives the neighbours which strongly
% influence i AND are in the cNodes set (C_i)
cNodeLogical = zeros(1, nrVar);
cNodeLogical(cNodes) = 1;
coarseInterp = bsxfun(@and, auxStrength, cNodeLogical);

% fNodeStrong: ith row gives neighbours which strongly influence i but are not 
% in the cNodes set (hence they are in the fNodes set) (D^s_i)
fNodeLogical = zeros(1, nrVar);
fNodeLogical(fNodes) = 1;
fNodeStrong = bsxfun(@and, auxStrength, fNodeLogical);

% weakNeighbours: ith row gives neighbours which do not strongly influence i
% (the left over neighbours) (D^w_i)
zeroDiagAh = Ah - spdiags(diag(Ah), 0, nrVar, nrVar);
weakNeighbours = xor(logical(zeroDiagAh), (coarseInterp | fNodeStrong));

% general interpolation rule:
% (I^h_{2h}e)_i = e_i 							if i in cNodes, 
% 				= sum_{j in C_i} w_{ij} e_j 	if i in fNodes
% here C_i is defined by coarseInterp(i, :)

if strcmp(setup.amg.interpMethod, 'classical')
	% classical interpolation (requires strong F-F connections to have a common C-point)

	% requires some attention (too much is calculated now, and divided by zero)
	% storage of NaN?
	
	% precompute denominator: denom(i) = aii + sum_{n in D^w_i} a_{in}, 
	% for i = 1:nrVar, (result is a column vector)
	denom = diag(Ah) + sum(Ah.*weakNeighbours,2);

	% denom = diag(Ah);
	% for i = 1:nrVar
	% 	for n = 1:nrVar
	% 		if weakNeighbours(i, n)
	% 			denom(i) = denom(i) + Ah(i,n);
	% 		end
	% 	end
	% end

	% precompute denominator in numerator: denomNum(i,m) = sum_{k in C_i} a_{mk},
	% for i,m = 1:nrVar (result is a matrix)
	denomNum = coarseInterp*Ah';

	% denomNum = zeros(nrVar);
	% for i = 1:nrVar
	% 	for m = 1:nrVar
	% 		for k = 1:nrVar
	% 			if coarseInterp(i, k)
	% 				denomNum(i,m) = denomNum(i,m) + Ah(m,k); 
	% 			end
	% 		end
	% 	end
	% end

	% compute numerator: 
	% numer(i,j) = a_{ij} + sum_{m in D^s_i (a_im amj)/denomNum(i,m)} 
	numer = Ah(:, :) + fNodeStrong.*(Ah./denomNum)*Ah(:, :);

	% numer = Ah;
	% for i = 1:nrVar
	% 	for j = 1:nrVar
	% 		for m = 1:nrVar
	% 			if fNodeStrong(i, m)
	% 				numer(i,j) = numer(i,j) + Ah(i, m)*Ah(m, j)/denomNum(i,m);
	% 			end
	% 		end
	% 	end
	% end

	% resulting weights: interpWeights(i,j) = - (numer(i,j)/denom(i))
	interpWeights = - bsxfun(@times, numer, 1./denom);

	% coarse node is interpolated by itself
	interpWeights(cNodes, cNodes) =  speye(length(cNodes));
	interpOp = interpWeights(:,cNodes);
	interpOp(isnan(interpOp)) = 0; % fix this!
elseif strcmp(setup.amg.interpMethod, 'F-F')
	% F-F interpolation 

end



end