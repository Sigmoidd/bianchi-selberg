import mpmath as mp
from screen import gfun, gpp0, field
mp.mp.dps = 30
# independent ground truth: int h r^2 dr  = 2*pi*(-g''(0))  -> compare with direct quadrature of sinc^{2k}
def I_direct(k, dl):
    f = lambda r: (mp.sin(dl*r)/(dl*r))**(2*k) * r*r
    # integrate far enough: tail ~ dl^{-2k} R^{3-2k}/(2k-3); use oscillatory quadrature via many panels
    return 2*mp.quad(f, mp.linspace(0, 400, 800)) / (2*mp.pi)   # = -g''(0) estimate (tail dropped)
def gpp_fd(k, dl):   # the repo's 4-point stencil
    hh = dl/2
    a0,a1,a2,a3 = (gfun(mp.mpf(0),k,dl), gfun(hh,k,dl), gfun(2*hh,k,dl), gfun(3*hh,k,dl))
    return (2*a0-5*a1+4*a2-a3)/hh**2
for k in (2,3,4):
    dl = mp.mpf('0.2')
    ex = -gpp0(k, dl); fd = -gpp_fd(k, dl)
    print(f"k={k}: exact -g''(0+)={mp.nstr(ex,12)}   repo 4-pt stencil={mp.nstr(fd,12)}   rel.err={mp.nstr((fd-ex)/ex,4)}")
# direct check of the exact value for k=2,3 (tail of order R^{3-2k} included only for k=3 roughly)
for k in (2,3):
    dl = mp.mpf('0.2')
    print(f"k={k}: direct quadrature (R=400) -g''(0) ~ {mp.nstr(I_direct(k,dl),10)}   exact formula {mp.nstr(-gpp0(k,dl),10)}")
