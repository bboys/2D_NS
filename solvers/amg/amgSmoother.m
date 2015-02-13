function u = amgSmoother(u, A, f, setup)
% one step of Gauss-Seidel
nrVar = size(A,1);
for row = nrVar
	u(row) = (f(row) - A(row,row + 1:end-1)*u(row + 1:end-1) -...
		A(row,2:row - 1)*u(2:row - 1))/A(row,row);
end


end