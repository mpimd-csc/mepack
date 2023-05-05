#!/bin/bash

#MISC
sed -f cmake/convert_double_single.sed src/misc/dla_small_solve8.f90 > src/misc/sla_small_solve8.f90
sed -f cmake/convert_double_single.sed src/misc/dla_small_solve4.f90 > src/misc/sla_small_solve4.f90
sed -f cmake/convert_double_single.sed src/misc/dla_dtrmm_fix.f90 > src/misc/sla_strmm_fix.f90
# sed -f cmake/convert_double_single.sed src/misc/dla_geaxb.f90 > src/misc/sla_geaxb.f90
sed -f cmake/convert_double_single.sed src/misc/dla_qtrmm2.f90 > src/misc/sla_qtrmm2.f90
sed -f cmake/convert_double_single.sed src/misc/dla_sort_ev.f90 > src/misc/sla_sort_ev.f90
sed -f cmake/convert_double_single.sed src/misc/dla_sort_gev.f90 > src/misc/sla_sort_gev.f90
sed -f cmake/convert_double_single.sed src/misc/dla_itranspose.f90 > src/misc/sla_itranspose.f90
sed -f cmake/convert_double_single.sed src/misc/dla_transform_standard.f90 > src/misc/sla_transform_standard.f90
sed -f cmake/convert_double_single.sed src/misc/dla_transform_general.f90 > src/misc/sla_transform_general.f90
sed -f cmake/convert_double_single.sed src/common/options/blocksize_double.f90 > src/common/options/blocksize_single.f90


# Recursive
sed -f cmake/convert_double_single.sed src/double/recursive/dla_trlyap_recursive.f90 > src/single/recursive/sla_trlyap_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_trsylv_recursive.f90 > src/single/recursive/sla_trsylv_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_trsylv2_recursive.f90 > src/single/recursive/sla_trsylv2_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_trstein_recursive.f90 > src/single/recursive/sla_trstein_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_tgsylv_recursive.f90 > src/single/recursive/sla_tgsylv_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_tglyap_recursive.f90 > src/single/recursive/sla_tglyap_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_tgstein_recursive.f90 > src/single/recursive/sla_tgstein_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_tgcsylv_recursive.f90 > src/single/recursive/sla_tgcsylv_recursive.f90
sed -f cmake/convert_double_single.sed src/double/recursive/dla_tgcsylv_dual_recursive.f90 > src/single/recursive/sla_tgcsylv_dual_recursive.f90


# Kernels
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv_kernel_44_nn.f90 > src/single/kernels/sla_trsylv_kernel_44_nn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv_kernel_44_nt.f90 > src/single/kernels/sla_trsylv_kernel_44_nt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv_kernel_44_tn.f90 > src/single/kernels/sla_trsylv_kernel_44_tn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv_kernel_44_tt.f90 > src/single/kernels/sla_trsylv_kernel_44_tt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv2_kernel_44_nn.f90 > src/single/kernels/sla_trsylv2_kernel_44_nn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv2_kernel_44_nt.f90 > src/single/kernels/sla_trsylv2_kernel_44_nt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv2_kernel_44_tn.f90 > src/single/kernels/sla_trsylv2_kernel_44_tn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trsylv2_kernel_44_tt.f90 > src/single/kernels/sla_trsylv2_kernel_44_tt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgsylv_kernel_44_nn.f90 > src/single/kernels/sla_tgsylv_kernel_44_nn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgsylv_kernel_44_nt.f90 > src/single/kernels/sla_tgsylv_kernel_44_nt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgsylv_kernel_44_tn.f90 > src/single/kernels/sla_tgsylv_kernel_44_tn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgsylv_kernel_44_tt.f90 > src/single/kernels/sla_tgsylv_kernel_44_tt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trlyap_kernel_44_n.f90 > src/single/kernels/sla_trlyap_kernel_44_n.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_trlyap_kernel_44_t.f90 > src/single/kernels/sla_trlyap_kernel_44_t.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgsylv_kron_kernel.f90 > src/single/kernels/sla_tgsylv_kron_kernel.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_kernel_44_nn.f90 > src/single/kernels/sla_tgcsylv_kernel_44_nn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_kernel_44_nt.f90 > src/single/kernels/sla_tgcsylv_kernel_44_nt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_kernel_44_tn.f90 > src/single/kernels/sla_tgcsylv_kernel_44_tn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_kernel_44_tt.f90 > src/single/kernels/sla_tgcsylv_kernel_44_tt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_dual_kernel_44_nn.f90 > src/single/kernels/sla_tgcsylv_dual_kernel_44_nn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_dual_kernel_44_nt.f90 > src/single/kernels/sla_tgcsylv_dual_kernel_44_nt.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_dual_kernel_44_tn.f90 > src/single/kernels/sla_tgcsylv_dual_kernel_44_tn.f90
sed -f cmake/convert_double_single.sed src/double/kernels/dla_tgcsylv_dual_kernel_44_tt.f90 > src/single/kernels/sla_tgcsylv_dual_kernel_44_tt.f90



#level 2
sed -f cmake/convert_double_single.sed src/double/level2/dla_trlyap_l2.f90 > src/single/level2/sla_trlyap_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trlyap_l2_opt.f90 > src/single/level2/sla_trlyap_l2_opt.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trstein_l2.f90 > src/single/level2/sla_trstein_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tglyap_l2.f90 > src/single/level2/sla_tglyap_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tglyap_l2_unopt.f90 > src/single/level2/sla_tglyap_l2_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgstein_l2.f90 > src/single/level2/sla_tgstein_l2.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2.f90 > src/single/level2/sla_trsylv_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_unopt.f90 > src/single/level2/sla_trsylv_l2_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy.f90 >     src/single/level2/sla_trsylv_l2_opt_local_copy.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy_32.f90 >  src/single/level2/sla_trsylv_l2_opt_local_copy_32.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy_64.f90 >  src/single/level2/sla_trsylv_l2_opt_local_copy_64.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy_96.f90 >  src/single/level2/sla_trsylv_l2_opt_local_copy_96.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy_128.f90 > src/single/level2/sla_trsylv_l2_opt_local_copy_128.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_reorder.f90 >     src/single/level2/sla_trsylv_l2_opt_reorder.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv_l2_opt_local_copy_xx.f90 >     src/single/level2/sla_trsylv_l2_opt_local_copy_xx.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2.f90 > src/single/level2/sla_trsylv2_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_unopt.f90 > src/single/level2/sla_trsylv2_l2_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy.f90 >     src/single/level2/sla_trsylv2_l2_opt_local_copy.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy_32.f90 >  src/single/level2/sla_trsylv2_l2_opt_local_copy_32.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy_64.f90 >  src/single/level2/sla_trsylv2_l2_opt_local_copy_64.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy_96.f90 >  src/single/level2/sla_trsylv2_l2_opt_local_copy_96.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy_128.f90 > src/single/level2/sla_trsylv2_l2_opt_local_copy_128.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_reorder.f90 >     src/single/level2/sla_trsylv2_l2_opt_reorder.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_trsylv2_l2_opt_local_copy_xx.f90 >     src/single/level2/sla_trsylv2_l2_opt_local_copy_xx.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_gardiner_laub.f90 > src/single/level2/sla_tgsylv_gardiner_laub.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2.f90 > src/single/level2/sla_tgsylv_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_rowwise.f90 > src/single/level2/sla_tgsylv_l2_rowwise.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_colwise.f90 > src/single/level2/sla_tgsylv_l2_colwise.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy.f90 >     src/single/level2/sla_tgsylv_l2_opt_local_copy.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy_32.f90 >  src/single/level2/sla_tgsylv_l2_opt_local_copy_32.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy_64.f90 >  src/single/level2/sla_tgsylv_l2_opt_local_copy_64.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy_96.f90 >  src/single/level2/sla_tgsylv_l2_opt_local_copy_96.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy_128.f90 > src/single/level2/sla_tgsylv_l2_opt_local_copy_128.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_reorder.f90 >     src/single/level2/sla_tgsylv_l2_reorder.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgsylv_l2_opt_local_copy_xx.f90 >     src/single/level2/sla_tgsylv_l2_opt_local_copy_xx.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2.f90 > src/single/level2/sla_tgcsylv_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_unopt.f90 > src/single/level2/sla_tgcsylv_l2_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy.f90 >     src/single/level2/sla_tgcsylv_l2_opt_local_copy.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy_32.f90 >  src/single/level2/sla_tgcsylv_l2_opt_local_copy_32.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy_64.f90 >  src/single/level2/sla_tgcsylv_l2_opt_local_copy_64.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy_96.f90 >  src/single/level2/sla_tgcsylv_l2_opt_local_copy_96.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy_128.f90 > src/single/level2/sla_tgcsylv_l2_opt_local_copy_128.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_reorder.f90 >     src/single/level2/sla_tgcsylv_l2_opt_reorder.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_l2_opt_local_copy_xx.f90 >     src/single/level2/sla_tgcsylv_l2_opt_local_copy_xx.f90

sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2.f90 > src/single/level2/sla_tgcsylv_dual_l2.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy.f90 >     src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy_32.f90 >  src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy_32.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy_64.f90 >  src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy_64.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy_96.f90 >  src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy_96.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy_128.f90 > src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy_128.f90
sed -f cmake/convert_double_single.sed src/double/level2/dla_tgcsylv_dual_l2_opt_local_copy_xx.f90 >     src/single/level2/sla_tgcsylv_dual_l2_opt_local_copy_xx.f90







#level3
sed -f cmake/convert_double_single.sed src/double/level3/dla_trlyap_l3.f90 > src/single/level3/sla_trlyap_l3.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trlyap_l3_2stage.f90 > src/single/level3/sla_trlyap_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv_l3_opt.f90 > src/single/level3/sla_trsylv_l3_opt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv_l3_unopt.f90 > src/single/level3/sla_trsylv_l3_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv_l3_2stage.f90 > src/single/level3/sla_trsylv_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv2_l3_opt.f90 > src/single/level3/sla_trsylv2_l3_opt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv2_l3_unopt.f90 > src/single/level3/sla_trsylv2_l3_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trsylv2_l3_2stage.f90 > src/single/level3/sla_trsylv2_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trstein_l3.f90 > src/single/level3/sla_trstein_l3.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_trstein_l3_2stage.f90 > src/single/level3/sla_trstein_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgsylv_l3_opt.f90 > src/single/level3/sla_tgsylv_l3_opt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgsylv_l3_colwise.f90 > src/single/level3/sla_tgsylv_l3_colwise.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgsylv_l3_rowwise.f90 > src/single/level3/sla_tgsylv_l3_rowwise.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgsylv_l3_2stage.f90 > src/single/level3/sla_tgsylv_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tglyap_l3.f90 > src/single/level3/sla_tglyap_l3.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tglyap_l3_2stage.f90 > src/single/level3/sla_tglyap_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgstein_l3.f90 > src/single/level3/sla_tgstein_l3.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgstein_l3_2stage.f90 > src/single/level3/sla_tgstein_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgcsylv_l3_opt.f90 > src/single/level3/sla_tgcsylv_l3_opt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgcsylv_l3_unopt.f90 > src/single/level3/sla_tgcsylv_l3_unopt.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgcsylv_l3_2stage.f90 > src/single/level3/sla_tgcsylv_l3_2stage.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgcsylv_dual_l3.f90 > src/single/level3/sla_tgcsylv_dual_l3.f90
sed -f cmake/convert_double_single.sed src/double/level3/dla_tgcsylv_dual_l3_2stage.f90 > src/single/level3/sla_tgcsylv_dual_l3_2stage.f90



# DAG
sed -f cmake/convert_double_single.sed src/double/dag/dla_trlyap_dag.f90 > src/single/dag/sla_trlyap_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_trsylv_dag.f90 > src/single/dag/sla_trsylv_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_trsylv2_dag.f90 > src/single/dag/sla_trsylv2_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_trstein_dag.f90 > src/single/dag/sla_trstein_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_tgsylv_dag.f90 > src/single/dag/sla_tgsylv_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_tglyap_dag.f90 > src/single/dag/sla_tglyap_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_tgstein_dag.f90 > src/single/dag/sla_tgstein_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_tgcsylv_dag.f90 > src/single/dag/sla_tgcsylv_dag.f90
sed -f cmake/convert_double_single.sed src/double/dag/dla_tgcsylv_dual_dag.f90 > src/single/dag/sla_tgcsylv_dual_dag.f90



# Frontends
for F in dla_ggcsylv dla_ggcsylv_dual dla_gelyap dla_gestein dla_gesylv dla_gesylv2 dla_gglyap dla_ggstein dla_ggsylv
do
	FS=`echo $F | sed -e 's/^d/s/'`
	sed -f cmake/convert_double_single.sed src/double/frontend/${F}.f90 > src/single/frontend/${FS}.f90
done

# Refine
for F in dla_ggcsylv dla_ggcsylv_dual dla_gelyap dla_gestein dla_gesylv dla_gesylv2 dla_gglyap dla_ggstein dla_ggsylv
do
	FS=`echo $F | sed -e 's/^d/s/'`
	sed -f cmake/convert_double_single.sed src/double/refine/${F}_refine.f90 > src/single/refine/${FS}_refine.f90
done

# Tests
for v in cgsylv ggsylv gesylv gesylv2 gelyap gestein gglyap ggstein
do
	sed -f cmake/convert_double_single.sed -f cmake/convert_double_single_c.sed test/arguments/double/check_arguments_$v.c > test/arguments/single/check_arguments_single_$v.c
done

