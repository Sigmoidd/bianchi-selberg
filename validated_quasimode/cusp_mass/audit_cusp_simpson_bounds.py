"""Exact rational audit of *global* Simpson and cusp truncation error bounds.

No special-function, floating-point or external symbolic package required.
The inequalities here are also listed explicitly in the proof note.
"""
from fractions import Fraction as F
from math import factorial
import json

PI_LO=F(314,100)
PI_HI=F(22,7)
ZMIN=2*PI_LO*F(5,3)
ZMAX=2*PI_HI*F(21,10)
assert ZMIN == F(157,15) and ZMAX < 14
assert F(6622119348586011,10**15) < 7

def exp_taylor(t,n):
    return sum((t**k/F(factorial(k))) for k in range(n+1))
assert exp_taylor(ZMIN,14)>30000
assert exp_taylor(F(3),8)>20
# e^(-z*cosh(3))/(z*sinh 3), cosh3 >10, sinh3 >9.
assert F(1,30000**10*90)<F(1,10**40)
D1=14
D2=14**2+14
D3=14**3+3*14**2+14
D4=14**4+6*14**3+7*14**2+14
P4=D4+4*7*D3+6*7**2*D2+4*7**3*D1+7**4
Mt=F(3,5)*F(P4,30000)
assert Mt<5, Mt
My=F(3,30000)*(7**4+4*7**3+12*7**2+24*7+24)
assert My<F(1,2),My
maxwidth=F(1,5)
ht=F(3,512)
hy=maxwidth/32
err_t=maxwidth*3*ht**4*5/180
err_y=maxwidth*hy**4*F(1,2)/180
err_tail=F(1,10**40)
err_total=err_t+err_y+err_tail
assert err_total<F(1,10**10)
report={'tier':'exact-rational proof audit',
'z_min':str(ZMIN),'z_max':str(ZMAX),
'exp_lower_14_greater_30000':True,'exp3_lower_8_greater_20':True,
'fourth_t_derivative_bound_upper':str(Mt),'fourth_y_derivative_bound_upper':str(My),
'error_t_simpson_under':str(float(err_t)),
'error_y_simpson_under':str(float(err_y)),
'error_truncated_t_under':'1e-40',
'error_total_less_than':'1e-10',
'checks_passed':8}
print(json.dumps(report,indent=2))
with open('cusp_error_audit.json','w') as fh:json.dump(report,fh,indent=2)
