function v = preconFunc(vold, setup, precon, Q,preconType);
% performs preconditioning step
v = vold;


if preconType == 2
	% ichol/ilu block preconditioner
	nrp = size(Q,1);
	v = precon\vold;
	v = precon'\v;
elseif preconType == 1
	% use AMG block preconditioner

	v(1:precon.nrv, :) = amgSolve([], vold(1:precon.nrv,:),...
		zeros(precon.nrv, size(vold,2)), setup, precon);

	v(precon.nrv + 1:end, :) = vold(precon.nrv + 1:end, :)./diag(Q);
	
	% idrstab (not working)
	% v(precon.nrv + 1:end, :) = bsxfun(@times, vold(precon.nrv + 1:end, :),...
	% 	1./diag(Q));

end
% if neither then do nothing (no precon)
end