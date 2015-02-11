maxit = 1500;
tol = 1e-8;

% test for s = S1:Ssn, l = L1:Lln
saveFile = 0;
sn = default('Number of s', 5);
ln = default('Number of l', 5);

S = 1 + [0:sn-1];
L = 1 + [0:ln-1];


lidVel = 1;
% create local matrices
[ localMatrix, basisOrder ] = createRectBasis();
nrPBasisF = size(localMatrix.pdivv.x,1);

% create mesh 
feMesh = createRectMesh(basisOrder);
stabC = createStabC(feMesh, localMatrix);


nrElts = feMesh.problemSize(1)*feMesh.problemSize(2);
nrNodes = feMesh.problemSize(3)*feMesh.problemSize(4);

% assemble matrices
% M = vmassAssembly( feMesh, localMatrix.vmass); % mass matrix
globalMatrix.L = PdivVAssembly( feMesh, localMatrix.pdivv); % "pressure mass" matrix
% D = diffusionAssembly( feMesh, localMatrix.stiff); % diffusive matrix
globalMatrix.D = laplaceAssembly( feMesh, localMatrix.stiff); % alternative (not using Sij)


% determine free nodes (interior)
% lidNodes = unique(feMesh.boundary.gamma2(:)); % requires regularisation
lidNodes = feMesh.boundary.gamma2(:, 2:end-1);
lidNodes = unique(lidNodes(:));

% homogeneous Dirichlet Nodes
homDNodes = unique([feMesh.boundary.gamma1(:);feMesh.boundary.gamma3(:);...
	feMesh.boundary.gamma4(:)]);

fixedVel = [lidNodes; lidNodes + nrNodes; homDNodes; homDNodes + nrNodes;];
freeVel = setdiff(1:2*nrNodes, fixedVel);

fixedPressure = [1];
freePressure = setdiff(1:nrPBasisF*nrElts, fixedPressure);

freeSol = [freeVel, 2*nrNodes + freePressure]; % include pressure DOF
fixedSol = [fixedVel', 2*nrNodes + fixedPressure];

nodeType = struct('freeVel', freeVel, 'freePressure', freePressure,...
	'freeSol', freeSol, 'fixedVel', fixedVel, 'fixedPressure', fixedPressure);

% fill in the boundary conditions
solVec = zeros(2*nrNodes + nrPBasisF*nrElts,1);
% solVec(lidNodes) = lidVel*(1 - feMesh.node(1, lidNodes).^4); % regularised
solVec(lidNodes) = lidVel;

% solve corresponding stokes problem as initial guess
M = [globalMatrix.D, -globalMatrix.L'; -globalMatrix.L -stabC];
rhsVec = -M(:, fixedVel)*solVec(fixedVel);

% % store matrix for testing solve methods
A = globalMatrix.D(freeVel, freeVel);
B = globalMatrix.L(freePressure, freeVel);
Q = pmassAssembly(feMesh, localMatrix.pmass);
Q = Q(freePressure, freePressure);
C = stabC(freePressure, freePressure);
f = rhsVec(freeSol);
sol = solVec(freeSol);



n = length(f);



tic
nrp = size(Q,1);
nrv = size(A,1);

if ~exist('C')
	C = sparse(nrp, nrp); % stabC
end

setup.type='nofill'; 
P1 = ichol(A, setup);
% P1 = sqrt(spdiags(diag(A),0,nrv,nrv));
P2 = sqrt(spdiags(diag(Q),0,nrp,nrp));
P = [P1 sparse(nrv, nrp); sparse(nrp, nrv) P2];
Mnotsym = [A, -B'; B, -C];
Msym = [A, -B'; -B, -C];
% showResult('build precon', toc, sol2, Msym, f, 0)

% tic
% sol = Msym\f;
% showResult('backslash', toc, sol2, Msym, f, 0)


tic
[sol2, ~, ~, iter, resminres] = minres(Msym,f, tol, maxit,P, P');
time.minres = toc;

time.idrns = zeros(sn,1);
time.idrs = zeros(sn,1);
time.idrstabns = zeros(sn,ln); 
time.idrstabs = zeros(sn,ln); 

MV.idrns = zeros(sn,1);
MV.idrs = zeros(sn,1);
MV.idrstabns = zeros(sn,ln); 
MV.idrstabs = zeros(sn,ln); 

relres.idrns = zeros(sn,1);
relres.idrs = zeros(sn,1);
relres.idrstabns = zeros(sn,ln); 
relres.idrstabs = zeros(sn,ln); 
relres.minres = resminres;

sizeString = 0;
for idxs = 1:sn
	s = S(idxs);
	% Rs = orth(randn(size(f,1),s));
	Rs = orth(randn(size(f,1),s));


	for idxl = 1:ln
		l = L(idxl);

		fprintf(repmat('\b',sizeString,1));
		string = sprintf('s = %2.0f, l = %2.0f',[s,l]);
		sizeString = size(string);
		fprintf(string)

		
		tic
		[sol2,res] = IDRstab_precon(Mnotsym,f,tol,maxit,s,l,Rs,P,P',[]);
		time.idrstabns(idxs,idxl) = toc;
		MV.idrstabns(idxs, idxl) = length(res);
		residrstabns(idxs,idxl).res = res;
		if res(end)/norm(f) > tol
			time.idrstabns(idxs, idxl) = Inf;
		end

		% showResult('idrstab not sym', toc, sol2, Mnotsym, f, length(residrstabns))

		tic
		[sol2,res] = IDRstab_precon(Msym,f,tol,maxit,s,l,Rs,P,P',[]);
		time.idrstabs(idxs,idxl) = toc;
		MV.idrstabs(idxs, idxl) = length(res);
		residrstabs(idxs,idxl).res = res;
		if res(end)/norm(f) > tol
			time.idrstabs(idxs, idxl) = Inf;
		end

	end

	tic
	[sol2,res] = Bi_IDRs(Mnotsym,f,tol,maxit,s,Rs,0,P,P',[]);
	time.idrns(idxs) = toc;
	MV.idrns(idxs) = length(res);
	residrns(idxs).res = res;
	if res(end)/norm(f) > tol
		time.idrns(idxs) = Inf;
	end

	tic
	[sol2,res] = Bi_IDRs(Msym,f,tol,maxit,s,Rs,0,P,P',[]);
	time.idrs(idxs) = toc;
	MV.idrs(idxs) = length(res);
	residrs(idxs).res = res;
	if res(end)/norm(f) > tol
		time.idrs(idxs) = Inf;
	end

end
fprintf(repmat('\b',sizeString,1));
fprintf('\n')

[~, minidrs] = min(time.idrs);
[~, minidrns] = min(time.idrns);

minidrstabs  = [0 0]; % s,l
[~,idx]=min(time.idrstabs(:));
[minidrstabs(1),minidrstabs(2)]=ind2sub(size(time.idrstabs),idx);

minidrstabns  = [0 0]; % s,l
[~,idx]=min(time.idrstabns(:));
[minidrstabns(1),minidrstabns(2)]=ind2sub(size(time.idrstabns),idx);

[~, MVminidrs] = min(MV.idrs);
[~, MVminidrns] = min(MV.idrns);

MVminidrstabs  = [0 0]; % s,l
[~,idx]=min(MV.idrstabs(:));
[MVminidrstabs(1),MVminidrstabs(2)]=ind2sub(size(MV.idrstabs),idx);

MVminidrstabns  = [0 0]; % s,l
[~,idx]=min(MV.idrstabns(:));
[MVminidrstabns(1),MVminidrstabns(2)]=ind2sub(size(MV.idrstabns),idx);



if saveFile
	save(dataName, 'time', 'MV', 'relres') % avoid saving big matrices
end


figure(1)
semilogy(relres.minres/norm(f), 'k-')
hold on
semilogy((residrs(minidrs).res)/norm(f), 'b:')
semilogy((residrstabs(minidrstabs(1), minidrstabs(2)).res)/norm(f), 'r:')

semilogy((residrns(minidrns).res)/norm(f), 'b-')
semilogy((residrstabns(minidrstabns(1), minidrstabns(2)).res)/norm(f), 'r-')

hold off

legend('minres',['idr(',num2str(S(minidrs)),') sym'],...
	['idrstab(',num2str(S(minidrstabs(1))),',',num2str(L(minidrstabs(2))),') sym'],...
	 ['idr(',num2str(S(minidrns)),') not sym'],...
	  ['idrstab(',num2str(S(minidrstabns(1))),',',num2str(L(minidrstabns(2))),') not sym'])


figure(4)
semilogy(relres.minres/norm(f), 'k-')
hold on
semilogy((residrs(MVminidrs).res)/norm(f), 'b:')
semilogy((residrstabs(MVminidrstabs(1), MVminidrstabs(2)).res)/norm(f), 'r:')

semilogy((residrns(MVminidrns).res)/norm(f), 'b-')
semilogy((residrstabns(MVminidrstabns(1), MVminidrstabns(2)).res)/norm(f), 'r-')

hold off



legend('minres',['idr(',num2str(S(MVminidrs)),') sym'],...
['idrstab(',num2str(S(MVminidrstabs(1))),',',num2str(L(MVminidrstabs(2))),') sym'],...
 ['idr(',num2str(S(MVminidrns)),') not sym'],...
  ['idrstab(',num2str(S(MVminidrstabns(1))),',',num2str(L(MVminidrstabns(2))),') not sym'])


  
cmap = 'gray';

figure(2)
subplot(2,2,1)
imagesc(L, S, time.idrstabns)
xlabel('l')
ylabel('s')
title('not sym time')
colorbar
colormap(cmap);
hold on
plot(L(minidrstabns(2)), S(minidrstabns(1)), 'rx', 'markersize', 10, 'linewidth', 3)
hold off


subplot(2,2,2)
imagesc(L, S, time.idrstabs)
xlabel('l')
ylabel('s')
title('sym time')
colorbar
colormap(cmap);
hold on
plot(L(minidrstabs(2)), S(minidrstabs(1)), 'rx', 'markersize', 10, 'linewidth', 3)
hold off

subplot(2,2,3)
imagesc(L, S, MV.idrstabns)
xlabel('l')
ylabel('s')
title('not sym MV')
colorbar
colormap(cmap);

subplot(2,2,4)
imagesc(L, S, MV.idrstabs)
xlabel('l')
ylabel('s')
title('sym MV')
colorbar
colormap(cmap);

figure(3)
imagesc(L, S, 100*(time.idrstabs-time.idrstabns)./time.idrstabns)
xlabel('l')
ylabel('s')
colormap(cmap);
colorbar
% end

% function [] = showResult(name, time, sol2, M, f, iter)
% fprintf(['%3.3f seconds, relres = %3.2e, MV = %4.0f by: ',name, '\n'], [time, norm(f - M*sol2)/norm(f), iter])
% end