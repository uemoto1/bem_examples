from scipy.special import jv, yv, jvp, yvp
from numpy import vectorize, linspace, zeros_like, exp, cos, pi, real, imag

def hv(v, z):
  return jv(v, z)+1.0j*yv(v, z)

def hvp(v, z):
  return jvp(v, z)+1.0j*yvp(v, z)

def coeff(k, r, m):
  return -(1.0j**m)*jv(m, k*r)/hv(m, k*r)

def scat(k, r, lmax=10):
  s = 0.0
  for m in xrange(-lmax, lmax+1):
    s += (4.0/k)*abs(get_coeff(k, r, m))**2
  return s

vscat = vectorize(scat)

def field(xs, t, k, r, lmax=10):
  f = zeros_like(xs, dtype=complex)

  for (i, x) in enumerate(xs):
    if (x <= r):
      f[i] = 0.0
    else:
      f[i] = exp(1.0j*k*x*cos(t))

  for m in xrange(-lmax, lmax+1):
    c = coeff(k, r, m)
    for (i, x) in enumerate(xs):
      if r < x:
        f[i] += c*hv(m, k*x)*exp(1.0j*m*t)


  return f


if __name__ == "__main__":
  k = 1.0
  r = 1.0
  x = linspace(0.0, 5.0)
  y = field(x, 0.0, k, r);
  for i in xrange(len(x)):
    print("%f\t%e\t%e" % (x[i], real(y)[i], imag(y)[i]))
