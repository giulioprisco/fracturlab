"""Deterministic expectation of the old Gemini script (alpha=0.4, eps=0.05, 2000 steps).
E phi_eps(B_s) = (2 pi (s^{2H}+eps^2))^{-1/2}; path point k sits at s_k = k/1999 (fbm package, n=1999),
the script's plotting grid is t = linspace(1e-4, 1, 2000). Also the eps -> 0 (no kernel) version."""
import numpy as np
from scipy.special import gamma
alpha=0.4; H=1-alpha; n=2000; a,b=1.0,2.0
t=np.linspace(1e-4,1,n); dt=1/n; f=a+b*t
s=np.arange(n)/(n-1)
an=a*t**alpha/gamma(alpha+1)+b*t**(alpha+1)/gamma(alpha+2)
def run(eps):
    if eps==0:
        dens=np.where(s>0,1/np.sqrt(2*np.pi)*np.maximum(s,1e-300)**(-H),0.0)
    else:
        dens=1/np.sqrt(2*np.pi*(s**(2*H)+eps**2))
    num=np.cumsum(f*dens)*dt
    R=an[-1]/num[-1]; sc=num*R
    out={'R':R}
    for x in (0.1,0.25,0.5,0.75):
        k=int(x*n); out[x]=(sc[k]-an[k])/an[k]*100
    return out
# continuous versions (no grid): forward integral int_0^x f(u) u^{-H}/sqrt(2pi) du vs I^a f(x)
from scipy.integrate import quad
def cont(eps):
    g=lambda u:(a+b*u)/np.sqrt(2*np.pi*(u**(2*H)+eps**2))
    I=lambda x:a*x**alpha/gamma(alpha+1)+b*x**(alpha+1)/gamma(alpha+2)
    F=lambda x:quad(g,0,x,limit=200)[0]
    R=I(1)/F(1); return {'R':R, **{x:(R*F(x)/I(x)-1)*100 for x in (0.1,0.25,0.5,0.75)}}
for e in (0.05,0):
    print('grid eps',e,run(e))
    print('cont eps',e,cont(e))
