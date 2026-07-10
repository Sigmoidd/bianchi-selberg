import bianchi_omega as M
import mpmath as mp
print("Robustness of B<1 over test functions (k, frac):")
worst=mp.mpf(0)
for k in [2,3]:
    for frac in [0.999,0.9,0.8,0.7]:
        Bw=M.compute_B('omega',k=k,frac=frac,R=80,verbose=False)
        Bi=M.compute_B('i',k=k,frac=frac,R=80,verbose=False)
        tag = "OK" if Bw<1 else "FAIL"
        worst=max(worst,Bw)
        print(f"  k={k} frac={frac}: B(Zi)={mp.nstr(Bi,6):>10}  B(Zomega)={mp.nstr(Bw,6):>10}  [{tag}]")
print(f"\nworst-case B(Zomega) = {mp.nstr(worst,6)}  ({'all < 1' if worst<1 else 'SOME >= 1'})")
