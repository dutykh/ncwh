# Matrices assembling for a classical Schwarzschild wormhole' quasi-normal modes computation
# Non-extreme case: BC1
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> 1/8*(1 + 2*s)*(s - 1)*(1 - y) - 1/2*L*(L+1):
  L01 := y -> -1/8*(1 - y)*(5*y + 3):
  L02 := y -> 1/4*(1 + y)*(1 - y)^2:
  L10 := y -> 3/16*(1 - y):
  L11 := y -> -1/4*(1 + y)*(y - 5):
  L12 := y -> 0:
  L20 := y -> 1/16*(7 - y):
  L21 := y -> 0:
  L22 := y -> 0:
  P := y -> simplify(add(a[j]*ChebyshevT(j, y), j=0..n-1)):
  M0 := Matrix(n):
  M1 := Matrix(n):
  M2 := Matrix(n):
  for i from 1 to n do
    xi := cos((2.0*i-1.0)*Pi/(2.0*n)); # Chebyshev roots collocation points
    expr0 := evalf(L00(xi)*P(xi) + L01(xi)*subs(x=xi, diff(P(x),x)) + L02(xi)*subs(x=xi, diff(P(x),x$2))):
    expr1 := evalf(L10(xi)*P(xi) + L11(xi)*subs(x=xi, diff(P(x),x)) + L12(xi)*subs(x=xi, diff(P(x),x$2))):
    expr2 := evalf(L20(xi)*P(xi) + L21(xi)*subs(x=xi, diff(P(x),x)) + L22(xi)*subs(x=xi, diff(P(x),x$2))):
    for j from 1 to n do
      M0[i,j] := coeff(expr0, a[j-1]):
      M1[i,j] := coeff(expr1, a[j-1]):
      M2[i,j] := coeff(expr2, a[j-1]):
    end do:
  end do:
  # We finally export the data from Maple and save in files:
  path := cat(p, "/data/"):
  nstr := convert(n, string);
  ExportMatrix(cat(path, "M0_", nstr, ".mat"), M0, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M1_", nstr, ".mat"), M1, target=MATLAB, mode=ascii):
  ExportMatrix(cat(path, "M2_", nstr, ".mat"), M2, target=MATLAB, mode=ascii):
end proc:
