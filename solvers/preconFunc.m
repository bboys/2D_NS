function v = preconFunc(vold, setup, precon, Q);
% performs preconditioning step
if strcmp(class(precon), 'double')
	% ichol block preconditioner
	nrp = size(Q,1);
	v = precon\vold;
	v = precon'\v;
elseif strcmp(class(precon), 'struct')
	v = vold;
	% use AMG block preconditioner
	v(1:precon.nrv) = amgSolve([], vold(1:precon.nrv),...
		zeros(precon.nrv,1), setup, precon);

	v(precon.nrv + 1:end) = diag(Q).\vold(precon.nrv + 1:end);
end

end