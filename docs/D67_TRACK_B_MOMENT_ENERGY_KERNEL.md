# Exact constant-kernel check

For strictly positive local B_T, zero finite local energy means each CR
interpolant is constant on its tetrahedron. Every local face mean is then
that constant. Shared identity rows equate the corresponding local constants.
An averaged row has positive weights summing to one; if all of its master
means already have the same constant, it equates its own tetrahedron with
that constant as well.

moment_kernel_check.cpp applies these exact implications on independently
verified moment rows. Starting with 9,697 identity-connected classes, two
passes merge all 909,276 tetrahedra into one class. Therefore the global
local-CR energy has only the constant kernel. This establishes no positive
spectral gap and does not certify the threshold matrix.
