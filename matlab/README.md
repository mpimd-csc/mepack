MATLAB / GNU Octave Interface                                         {#matlab}
=============================

MEPACK provides an interface for MATLAB and GNU Octave. The interface is
optional and its build need to be enabled separately. See
[Installation](../doc/install.md) for details. Once upon the interface is build,
the files either reside in
```shell
 ${BUILDDIR}/matlab/matlab
```
for the MATLAB interface and in
```shell
 ${BUILDDIR}/matlab/octave
```
for the GNU Octave interface. These paths need to be added to your MATLAB or
GNU Octave search path.

Interface Functions
-------------------

The MATLAB and GNU Octave interface provides access to all basic routines of
MEPACK. All function provide a detailed help page, which is accessible via
```matlab
  help mepack_XXXXXX
```

The following interfaces exist for general coefficient matrices:

 - `mepack_csylv`: Solve the generalized coupled Sylvester equation (`CSYLV`)
                   $`(AR \pm LB=E,CR \pm LD = F)`$.
 - `mepack_csylv_dual`: Solve the generalized coupled Sylvester equation
                        (`CSYLV_DUAL`) $`(AR + CL = E, \pm RB \pm LD =F)`$.
 - `mepack_glyap`: Solve the generalized Lyapunov equation (`GLYAP`)
                   $`AXB^{T} + BXA^T = Y`$.
 - `mepack_gstein`: Solve the generalized Stein equation (`GSTEIN`)
                    $`AXA^T - EXE^T = Y`$.
 - `mepack_gsylv`: Solve the generalized Sylvester equation (`GSYLV`)
                   $`AXB \pm CXD = Y`$.
 - `mepack_lyap`: Solve the standard Lyapunov equation (`LYAP`)
                  $`AX + XA^T = Y`$.
 - `mepack_stein`: Solve the standard Stein equation (`STEIN`)
                   $`AXA^T - X = Y`$.
 - `mepack_sylv2`: Solve the discrete-time Sylvester equation (`SYLV2`)
                   $`AXB - X = Y`$.
 - `mepack_sylv`: Solve the standard Sylvester equation (`SYLV`)
                  $`AX\pm XB = Y`$.

The following interfaces are available for (quasi-) triangular coefficient
matrices:

 - `mepack_tgcsylv`: Solve the generalized coupled Sylvester equation (`CSYLV`)
                     $`(AR \pm LB=E,CR \pm LD = F)`$.
 - `mepack_tgcsylv_dual`: Solve the generalized coupled Sylvester equation
                          (`CSYLV_DUAL`) $`(AR + CL = E, \pm RB \pm LD =F)`$.
 - `mepack_tglyap`: Solve the generalized Lyapunov equation (`GLYAP`)
                    $`AXB^{T} + BXA^T = Y`$.
 - `mepack_tgstein`: Solve the generalized Stein equation (`GSTEIN`)
                     $`AXA^T - EXE^T = Y`$.
 - `mepack_tgsylv`: Solve the generalized Sylvester equation (`GSYLV`)
                    $`AXB \pm CXD = Y`$.
 - `mepack_trlyap`: Solve the standard Lyapunov equation (`LYAP`)
                    $`AX + XA^T = Y`$.
 - `mepack_trstein`: Solve the standard Stein equation (`STEIN`)
                     $`AXA^T - X = Y`$.
 - `mepack_trsylv2`: Solve the discrete-time Sylvester equation (`SYLV2`)
                     $`AXB - X = Y`$.
 - `mepack_trsylv`: Solve the standard Sylvester equation (`SYLV`)
                    $`AX\pm XB = Y`$.

The following routines implement the iterative refinement strategy:

 - `mepack_csylv_refine`:  Solve the generalized coupled Sylvester equation
                           (`CSYLV`) $`(AR \pm LB=E,CR \pm LD = F)`$.
 - `mepack_csylv_dual_refine`: Solve the generalized coupled Sylvester equation
                            (`CSYLV_DUAL`) $`(AR + CL = E, \pm RB \pm LD =F)`$.
 - `mepack_glyap_refine`: Solve the generalized Lyapunov equation (`GLYAP`)
                          $`AXB^{T} + BXA^T = Y`$.
 - `mepack_gstein_refine`: Solve the generalized Stein equation (`GSTEIN`)
                           $`AXA^T - EXE^T = Y`$.
 - `mepack_gsylv_refine`: Solve the generalized Sylvester equation (`GSYLV`)
                          $`AXB \pm CXD = Y`$.
 - `mepack_lyap_refine`:  Solve the standard Lyapunov equation (`LYAP`)
                          $`AX + XA^T = Y`$.
 - `mepack_stein_refine`: Solve the standard Stein equation (`STEIN`)
                          $`AXA^T - X = Y`$.
 - `mepack_sylv2_refine`: Solve the discrete-time Sylvester equation (`SYLV2`)
                          $`AXB - X = Y`$.
 - `mepack_sylv_refine`: Solve the standard Sylvester equation (`SYLV`)
                         $`AX\pm XB = Y`$.


Known Issues
------------

### MATLAB and OpenMP Codes

Due to bugs in MATLABs GNU OpenMP - Intel OpenMP wrapper, algorithms using
OpenMP's task dependencies are disable when calling the interface from MATLAB.
The algorithms can be enabled at runtime by setting the `openmp` value in the
optional **options** structure to `1`.


