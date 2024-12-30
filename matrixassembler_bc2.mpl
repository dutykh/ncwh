# Matrices assembling for non-commutative Schwarzschild wormhole' quasi-normal modes computation
# Non-extreme case 2
MatrixAssembler := proc (
  d::integer, # d : number of digits used in computations
  n::integer, # n : number of Tchebyshev modes
  s::numeric, # s : spin (0, 1, 2)
  L::numeric, # L : angular momentum, L >= s
  mu::numeric, # mu : dimensionless BH's mass
  p::string   # p : string containing the path where we save the assembled matrices
  )
  if L < s then
    error(1, "Angular momentum cannot be smaller than spin!");
  end if:
  local f::function, g::function, xh::numeric, alpha::numeric, fp::function, gp::function, V::function, F::function, Fp::function, P::function, L00::function, L01::function, L02::function, L10::function, L11::function, L12::function, L20::function, L21::function, L22::function, M0::Matrix, M1::Matrix, M2::Matrix, i::integer, j::integer, expr0::algebraic, expr1::algebraic, expr2::algebraic, xi::numeric, path::string, nstr::string:
  with(LinearAlgebra):
  Digits := d:
  # Determination of the event horizon:
  f := x -> 1 - erf(mu*x)/x + 2*mu/sqrt(Pi)*exp(-mu^2*x^2):
  x0 := fsolve(f(x) = 0, x = 0.99):
  printf("Found WH's external throat x0 = %f\n", x0);
  # Other parameters appearing in the model:
  f := y -> 1 - 1/(2*x0)*(1 - y)*erf(2*mu*x0/(1 - y)) + 2*mu*exp(-4*mu^2*x0^2/(1 - y)^2)/sqrt(Pi):
  fp := y -> erf(2*mu*x0/(1 - y))/(2*x0) - 2*exp(-4*mu^2*x0^2/(1 - y)^2)*mu/((1 - y)*sqrt(Pi)) - 16*mu^3*x0^2*exp(-4*mu^2*x0^2/(1 - y)^2)/((1 - y)^3*sqrt(Pi)):
  V := y -> (1 - y)^2/4*(L*(L + 1) + 1/2*(1 - s)*(1 + 2*s)*(1 - y)*fp(y)):
  # Definition of the 2nd order ODE coefficients:
  L00 := y -> (1 - y)^2*((1 + y)/4*fp(y) - 1/4*(5 + 3*y)*f(y)/(1 - y) - 4*(1 + y)^2*V(y)/(1 - y)^4)/(1 + y)^2:
  L01 := y -> (1 - y)*((-y^2 + 1)/2*fp(y) - (3*y + 1)*f(y))/(1 + y):
  L02 := y -> (1 - y)^2*f(y):
  L10 := y -> ((1 + y)*(4*x0 + 1 - y)*fp(y) + 4*(2*x0 - y)*f(y))/(4*(1 + y)):
  L11 := y -> (4*x0 + 1 - y)*f(y):
  L12 := y -> 0:
  L20 := y -> -1/4*f(y)*(4*x0/(1 - y) + 1)^2 + 4*x0^2/(1 - y)^2:
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