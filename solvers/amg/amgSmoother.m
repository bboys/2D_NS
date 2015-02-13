function u = amgSmoother(u, A, f, setup)
if setup.amg.smoothType == 1
	% one step of Gauss-Seidel
	nrVar = size(A,1);
	for row = 1:nrVar
		u(row) = (f(row) - A(row,row + 1:end)*u(row + 1:end) -...
			A(row,1:row-1)*u(1:row-1 ))/A(row,row);
	end
end

end