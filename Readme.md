Todo
- MEX: sparse MV in Fortran
- MEX: minres Fortran code schrijven (preconditioned)
- stokes solver
-- preconditioned MINRES
-- multigrid (algebraisch / geometrisch)
- triangulation
- error estimation
- time stepping
- inexact newton solver
- quadrature
- amg
-- residual check toevoegen (zonder extra MV)

Problemen
- AMG interpolatie operator verbeteren




Verbeteringen
- vervang setdiff in createInterpOp met logical indexing
- nonlin1 is voor 2e term, nonlin2 voor 1e term? omdraaien?
- newton: als residu groter wordt, oplossing niet updaten als newton gebruikt wordt
- other nodal points
