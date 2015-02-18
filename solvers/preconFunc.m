function [v] = preconFuncStokes(vold, setup, precon);
% performs preconditioning step
v = vold;

if ~isempty(precon)
	if precon.type == 2
		% ichol/ilu block preconditioner
		nrp = size(precon.Q,1);
		v = precon.L\vold;
		v = precon.U\v;
	elseif precon.type == 1
		% use AMG block preconditioner
		v(1:precon.nrv, :) = amgSolve([], vold(1:precon.nrv,:),...
			zeros(precon.nrv, size(vold,2)), setup, precon);

		v(precon.nrv + 1:end, :) = vold(precon.nrv + 1:end, :)./diag(precon.Q);
	end
end

end