#Use mpif90 located at /usr/local/packages/mvapich/1.1/intel-10.1/bin/mpif90
#i.e. mvapich 1.1 compiled with Intel compiler version 10.1 
#by adding following to .soft file:
#+intel-fc-10.1
#+mvapich-1.1-intel-10.1
#Otherwise THE CODE WILL RUN BUT STARS WILL DANCE IN  NON-PHYSIXIE WAY!


.SUFFIXES: .F90 .f .c .o

F90FILES= main.F90 advect_frac.F90 advect_bc.F90 slope.F90 flux.F90 flux_frac.F90 \
 source.F90 source_bc.F90 make_source_temp.F90 make_source_temp_bc.F90 vel.F90    \
 vel_bc.F90 stress.F90 visc.F90 visc_time.F90 comm.F90 delta.F90 comm_dir.F90     \
 drag.F90 drag_bc.F90 save_state.F90 state.F90 potential_solver.F90 bessel.F90    \
 helmadi.F90 swap_r_z.F90 swap_z_r.F90 tridagr.F90 tridagz.F90 realft.F90         \
 comm_r_sweep.F90 comm_z_sweep.F90 swap_phi_r.F90 swap_r_phi.F90 output_files.F90 \
 output_kernel_pio.F90 dodiag.1.F90 diag_find_rhomax.F90 diag_make_mask.F90       \
 diag_find_lpoints.F90 diag_momments.F90 diag_spins.F90 diag_gwave.F90            \
 diag_write.1.F90 dsummer.F90 potsetup.F90 tm.F90 sm.F90 bm.F90 elle.F90 ellf.F90 \
 gammln.F90 rd.F90 rf.F90 setup_frac.F90 initialize_frac.F90 ritecont_frac.F90    \
 scfin_frac.F90 debugin.F90 calc_com_coords.F90 com_vel_accn.F90 ritecont_frac_new.F90

OFILES= $(F90FILES:.F90=.o) 

hydro:$(OFILES)

#	cmpif90c -O3 -o hydro $(OFILES)
#       mpif90 -O3 -o hydro -v $(OFILES)
	mpif90 -O0 -o hydro $(OFILES)

$(OFILES): runhydro.h

.F90.o: runhydro.h

# normal compilation after SP networkupgrade, don't use -qhot and -qipa
#	mpif90 -c -O3 -r8 -v $<
#	cmpif90c -c -O3 -r8 $<
#        mpif90 -c -O3 -r8 -v $<
	mpif90 -c -O0 -r8 $<

clean:
	/bin/rm -f *.o *.il hydro F*.f

