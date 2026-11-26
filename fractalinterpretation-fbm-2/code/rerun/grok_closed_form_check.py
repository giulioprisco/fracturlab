"""v0.4 check of the Grok script's printed 'proportional fractional integral' (H=0.15, alpha=0.85, f=1+t, T=1).
Shows that its closed form equals the left-sided RL integral I^alpha f(T), while 'analytical' is the
right-sided functional at 0, (2pi)^{-1/2} int_0^T f(u) u^{alpha-1} du = Gamma(alpha)/sqrt(2pi) * I^alpha_{T-} f(0)."""
from math import gamma, sqrt, pi
from scipy.integrate import quad
H = 0.15; al = 1 - H; a = b = 1.0; T = 1.0
code = (a + b*T)*T**al/(gamma(al)*al) - b*T**(al+1)/(gamma(al)*(al+1))
left_quad = quad(lambda s: (T-s)**(al-1)*(a+b*s), 0, T, limit=200)[0]/gamma(al)
left_closed = a*T**al/gamma(al+1) + b*T**(al+1)/gamma(al+2)
right_quad = quad(lambda u: u**(al-1)*(a+b*u), 0, T, limit=200)[0]/gamma(al)
print("script closed form            ", code)
print("left-sided I^a f(T), quadrature", left_quad)
print("left-sided I^a f(T), closed    ", left_closed)
print("right-sided functional at 0    ", right_quad)
print("Gamma/sqrt(2pi) * left  =", gamma(al)/sqrt(2*pi)*code, "(script prints 0.7230)")
print("Gamma/sqrt(2pi) * right =", gamma(al)/sqrt(2*pi)*right_quad, "(script 'analytical' 0.6850)")
