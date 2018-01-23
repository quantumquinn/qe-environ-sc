#!/bin/bash

# plugin_int_forces

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
  USE environ_base,  ONLY : env_static_permittivity, env_dielectric_regions, rhopol \
  USE environ_base,  ONLY : env_external_charges, rhoexternal \
  USE environ_main,  ONLY : calc_fenviron \
!Environ patch
' plugin_int_forces.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP), ALLOCATABLE :: force_environ(:,:), force_tmp(:,:) \
!Environ patch
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF (use_environ) THEN \
    ! \
    ALLOCATE(force_tmp(3,nat)) \
    ALLOCATE(force_environ(3,nat)) \
    ! \
    force_environ=0.0_dp \
    ! \
    IF(do_comp_mt) THEN \
      force_tmp=0.0_dp \
      ALLOCATE(auxr(dfftp%nnr)) \
      ALLOCATE(auxg(ngm)) \
      auxg = CMPLX(0.0_dp,0.0_dp) \
      auxr = CMPLX(0.0_dp,0.0_dp) \
      if(env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0) & \
        auxr = CMPLX(rhopol(:),0.0, kind=DP) \
      if(env_external_charges .GT. 0) & \
        auxr = auxr + CMPLX(rhoexternal(:),0.0, kind=DP) \
      CALL fwfft ("Dense", auxr, dfftp) \
      auxg(:)=auxr(nl(:)) \
      CALL wg_corr_force(.false.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, & \
                        1, auxg, force_tmp) \
      force_environ = force_environ + force_tmp \
      DEALLOCATE(auxr,auxg) \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9001 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! ... Computes here the solvent contribution \
    ! \
    IF ( env_static_permittivity .GT. 1.D0 .OR. env_dielectric_regions .GT. 0 ) THEN \
      force_tmp=0.0_dp \
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, & \
                     g, rhopol, nl, 1, gstart, gamma_only, vloc, & \
                     force_tmp ) \
      force_environ = force_environ + force_tmp \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9002 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! \
    ! ... Computes here the external charges contribution \
    ! \
    IF ( env_external_charges .GT. 0 ) THEN \
      force_tmp = 0.0_DP \
      CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, & \
                   g, rhoexternal, nl, 1, gstart, gamma_only, vloc, & \
                   force_tmp ) \
      force_environ = force_environ + force_tmp \
      IF ( iverbosity > 0 ) THEN \
        WRITE( stdout, 9003 ) \
        DO na = 1, nat \
           WRITE( stdout, 9035 ) na, ityp(na), ( force_tmp(ipol,na), ipol = 1, 3 ) \
        END DO \
      ENDIF \
    ENDIF \
    ! \
    ! ... Add the other environment contributions \
    ! \
    CALL calc_fenviron( dfftp%nnr, nspin, nat, force_environ ) \
    ! \
    DEALLOCATE(force_tmp) \
    ! \
    IF ( iverbosity > 0 ) THEN \
      WRITE( stdout, 9004 ) \
      DO na = 1, nat \
         WRITE( stdout, 9035 ) na, ityp(na), ( force_environ(ipol,na), ipol = 1, 3 ) \
      END DO \
      WRITE( stdout, * ) \
    ENDIF \
    ! \
    force = force_environ \
    ! \
    DEALLOCATE(force_environ) \
  END IF \
  ! \
9001 FORMAT(5x,"The MT-Environ correction contrib. to forces") \
9002 FORMAT(5x,"The dielectric solvent contribution to forces") \
9003 FORMAT(5x,"The external charges contribution to forces") \
9004 FORMAT(5x,"The global environment contribution to forces") \
9035 FORMAT(5X,"atom ",I4," type ",I2,"   force = ",3F14.8) \
!Environ patch
' tmp.2 > tmp.1

mv tmp.1 plugin_int_forces.f90

# plugin_read_input

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_input, ONLY : read_environ \
!Environ patch
' plugin_read_input.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL read_environ(nat, ntyp, assume_isolated, ibrav) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_read_input.f90

# plugin_clean

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_init, ONLY : environ_clean \
!Environ patch
' plugin_clean.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clean(lflag) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clean.f90

# plugin_summary

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_info, ONLY : environ_summary \
!Environ patch
' plugin_summary.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_summary() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_summary.f90

# plugin_initbase

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_base, ONLY : ir_end \
USE    environ_init, ONLY : environ_initbase \
!Environ patch
' plugin_initbase.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1)) \
  IF ( use_environ ) CALL environ_initbase( dfftp%nnr ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_initbase.f90

# plugin_clock

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_info, ONLY : environ_clock \
!Environ patch
' plugin_clock.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_clock(stdout) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_clock.f90

# plugin_print_energies

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE    environ_info, ONLY : environ_print_energies \
!Environ patch
' plugin_print_energies.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
   if(use_environ) CALL environ_print_energies() \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_print_energies.f90

# plugin_init_ions

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,            ONLY : alat, at \
USE ions_base,            ONLY : zv, nat, nsp, ityp, tau \
USE environ_init,         ONLY : environ_initions \
!Environ patch
' plugin_init_ions.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initions( dfftp%nnr, nat, nsp, ityp, zv, tau, alat, & \
                                             & at, tot_charge ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_ions.f90

# plugin_init_cell

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE cell_base,            ONLY : at, alat, omega, ibrav \
USE environ_init,         ONLY : environ_initcell \
!Environ patch
' plugin_init_cell.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF ( use_environ ) call environ_initcell( dfftp%nnr, dfftp%nr1, dfftp%nr2, dfftp%nr3, & \
                           ibrav, omega, alat, at ) \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_cell.f90

# plugin_scf_energy

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : deenviron, esolvent,    & \
                                 ecavity, epressure, eperiodic, eioncc,   & \
                                 eextcharge \
USE environ_main,          ONLY : calc_eenviron \
!Environ patch
' plugin_scf_energy.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) THEN \
       call calc_eenviron( dfftp%nnr, nspin, rhoin%of_r, deenviron, esolvent, & \
                            ecavity, epressure, eperiodic, eioncc, eextcharge ) \
        ! \
        plugin_etot = plugin_etot + deenviron + esolvent + ecavity + epressure + eperiodic + eioncc + eextcharge \
        ! \
  END IF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_energy.f90

# plugin_init_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : vltot_zero \
!Environ patch
' plugin_init_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  IF(use_environ) vltot_zero = vltot \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_init_potential.f90

# plugin_scf_potential

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_base,         ONLY : update_venviron, vltot_zero,    & \
                                 environ_thr, environ_restart \
USE environ_main,          ONLY : calc_venviron \
!Environ patch
' plugin_scf_potential.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
        ! \
        ! environ contribution to the local potential \
        ! \
        vltot = vltot_zero \
        !vltot_zero = 0.0_dp \
        ! \
        IF ( dr2 .GT. 0.0_dp ) THEN \
           update_venviron = .NOT. conv_elec .AND. dr2 .LT. environ_thr \
        ! \
        ELSE \
           update_venviron = environ_restart \
           ! for subsequent steps of optimization or dynamics, compute \
           ! environ contribution during initialization \
           IF ( .NOT. environ_restart ) environ_restart = .TRUE. \
        ENDIF \
        ! \
        IF ( update_venviron ) WRITE( stdout, 9200 ) \
        CALL calc_venviron( update_venviron, dfftp%nnr, nspin, dr2, rhoin%of_r, vltot ) \
        ! \
        !CALL sum_vrs( dfftp%nnr, nspin, vltot, vrs, vrs) \
        ! \
9200 FORMAT(/"     add environment contribution to local potential") \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 plugin_scf_potential.f90

# makov_payne

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
USE environ_mp,  ONLY : environ_makov_payne \
!Environ patch
' makov_payne.f90 > tmp.1

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
     IF(use_environ) THEN \
       CALL environ_makov_payne( dfftp%nnr, nspin, rho%of_r, x0, etot ) \
     ENDIF \
!Environ patch
' tmp.1 > tmp.2

mv tmp.2 makov_payne.f90

#plugin_initialization

sed '/Environ MODULES BEGIN/ a\
!Environ patch \
  USE environ_base,  ONLY : semiconductor_region, electrode_charge, flatband_fermi \
  USE environ_base,  ONLY : flatband_calc \
  USE semiconductor, ONLY : insert_helmholtz, read_flatband_pot \
!Environ patch 
' plugin_initialization.f90 > tmp.1


sed '/Environ CALLS BEGIN/ a\
!Environ patch \
  ! \
 \
!***************************************************************************** \
! \
! This checks on whether semiconductor optimization is used and either starts \
! the initial calculation of flatband potential or reads flatband potential from \
! file according to user input \
! \
!***************************************************************************** \
 \
  IF (use_environ) THEN \
 \
     IF (semiconductor_region) THEN \
        CALL start_clock( "semiconductor" ) \
        nstep = 100 \
        SELECT CASE( TRIM(flatband_calc)) \
 \
        CASE("from_input") \
            tot_charge =0.7*electrode_charge \
            insert_helmholtz = .TRUE. \
            CALL read_flatband_pot() \
            WRITE( stdout, 1001)flatband_fermi,tot_charge \
        CASE("from_scratch") \
            tot_charge = 0.0 \
            insert_helmholtz = .FALSE. \
            WRITE( stdout, 1002) tot_charge \
        CASE DEFAULT \
 \
            call errore("semiconductor","flatband_calc = " // & \
                        & trim(flatband_calc) // " not_recognized",4) \
 \
        END SELECT \
        CALL stop_clock( "semiconductor" ) \
 \
     END IF \
 \
  END IF \
 \
1001 FORMAT(5x,//"*******************************************"//,& \
               &"     Flatband potential of ",F14.8," read from input."// & \
            &   "     Using initial charge of:  ",F14.8,//& \
              & "*******************************************") \
1002 FORMAT(5x,//"*******************************************"//,& \
                &"     Running initial calculation for flatband."//& \
             &   "     Using charge of: ",F14.8,//& \
                &"*******************************************") \
 \
!Environ patch 
' tmp.1 > tmp.2

mv tmp.2 plugin_initialization.f90

#plugin_semiconductor_update
sed '/Environ MODULES BEGIN/ a\
!Environ patch \
  USE environ_base,  ONLY : semiconductor_region, flatband_calc, flatband_fermi \
  USE environ_base,  ONLY : electrode_charge, chg_thr \
  USE environ_base,  ONLY : flatband_calc \
  USE environ_base,  ONLY : environ_unit \
  USE semiconductor, ONLY : insert_helmholtz, save_flatband_pot \
  USE semiconductor, ONLY : save_current_pot, bulk_potential \
  USE semiconductor, ONLY : area \
!Environ patch 
' plugin_semiconductor_update.f90 > tmp.1

sed '/Environ VARIABLES BEGIN/ a\
!Environ patch \
  REAL(DP)                  ::   cur_chg \
  REAL(DP)                  ::   prev_chg, prev_chg2 \
  REAL(DP)                  ::   cur_dchg \
  REAL(DP)                  ::   cur_fermi \
  REAL(DP)                  ::   prev_dchg \
  REAL(DP)                  ::   gamma_mult \
  REAL(DP)                  ::   prev_step_size \
  REAL(DP)                  ::   ss_chg, charge \
  INTEGER                   ::   chg_step, na \
  REAL(DP)                  ::   surf_area \
  REAL(DP)                  :: sqbohrtosqcm \
  REAL(DP)                  :: chg_per_area \
  REAL(DP)                  :: ss_chg_per_area \
  REAL(DP)                  :: ss_potential, total_potential \
  REAL(DP)                  :: dft_chg_max, dft_chg_min \
  REAL(DP)                  :: change_vec \
  REAL(DP)                  :: v_cut \
  REAL(DP)                  :: ionic_charge \
  LOGICAL                   :: converge \
!Environ patch 
' tmp.1 > tmp.2

sed '/Environ CALLS BEGIN/ a\
!Environ patch \
 \
!************************************************* \
! \
! This section designed to run after a call to electrons. Basically, it checks \
! whether the semiconductor charge has converged and then updates the relevant \
! quantities (tot_charge) accordingly \
! \
!************************************************* \
 \
  gamma_mult = 0.15 \
 \
 \
  converge = .TRUE. \
  ionic_charge = 0._DP \
  DO na = 1, nat \
     ionic_charge = ionic_charge + zv( ityp(na) ) \
  END DO \
 \
 \
 \
  IF (use_environ .AND. semiconductor_region) THEN \
    CALL start_clock( "semiconductor" ) \
    IF (TRIM(flatband_calc) .EQ. "from_scratch") THEN \
         chg_step = i_step - 1 \
    ELSE \
         chg_step = i_step \
    END IF \
    !! Initializing the constraints of possible DFT charges \
    ! Should probably be initialized at chg_step =1 but that seemed to be \
    ! creating some trouble possibly \
    IF (chg_step == 1) THEN \
       IF (electrode_charge > 0.0) THEN \
          dft_chg_max = 2.0*electrode_charge \
          dft_chg_min = 0.0 \
       ELSE \
          dft_chg_min = 2.0*electrode_charge \
          dft_chg_max = 0.0 \
       END IF \
 \
    END IF \
 \
 \
    IF (.NOT. insert_helmholtz) THEN \
        insert_helmholtz = .TRUE. \
        tot_charge = 0.7*electrode_charge \
        flatband_fermi = ef*rytoev \
        conv_ions = .FALSE. \
        CALL qexsd_set_status(255) \
        CALL punch( "config" ) \
        CALL add_qexsd_step(i_step) \
        CALL save_flatband_pot(dfftp%nnr) \
        WRITE( stdout, 1001) flatband_fermi,tot_charge \
        ! \
        ! ... re-initialize atomic position-dependent quantities \
        ! \
        nelec = ionic_charge - tot_charge \
        CALL update_pot() \
        CALL hinit1() \
    ELSE \
         cur_fermi = ef*rytoev \
         CALL save_current_pot(dfftp%nnr,cur_fermi,cur_dchg,ss_chg,v_cut,chg_step) \
         !IF (ionode) THEN \
         ! making sure constraints are updated \
         IF (electrode_charge > 0) THEN \
            IF (v_cut >= 0.0 ) THEN \
                dft_chg_max = tot_charge \
                converge = .FALSE. \
            ELSE IF (ss_chg < 0.0) THEN \
                dft_chg_min = tot_charge \
                converge = .FALSE. \
            ELSE \
                prev_chg2 = tot_charge \
            END IF \
         ELSE \
            IF (v_cut <= 0.0 ) THEN \
                dft_chg_min = tot_charge \
                converge = .FALSE. \
            ELSE IF (ss_chg > 0.0) THEN \
                dft_chg_max = tot_charge \
                converge = .FALSE. \
            ELSE \
                prev_chg2 = tot_charge \
            END IF \
         END IF \
         CALL mp_bcast(dft_chg_min, ionode_id,intra_image_comm) \
         CALL mp_bcast(dft_chg_max, ionode_id,intra_image_comm) \
         IF (chg_step > 1 )THEN \
             gamma_mult = (cur_chg - prev_chg)/(cur_dchg - prev_dchg) \
         END IF \
         WRITE(environ_unit,*)"cur_chg: ",cur_chg \
         WRITE(environ_unit,*)"prev_chg: ",prev_chg \
         WRITE(environ_unit,*)"cur_dchg: ",cur_dchg \
         WRITE(environ_unit,*)"prev_dchg: ",prev_dchg \
         WRITE(environ_unit,*)"Using gamma of ",gamma_mult \
         change_vec = -gamma_mult*cur_dchg \
         prev_chg = tot_charge \
         ! This is my way of trying to impose limited constraints with an \
         ! unknown constraining function. Theres almost certainly a more \
         ! efficient way to do this but I havent thought of it yet \
 \
         IF ((tot_charge + change_vec) > dft_chg_max ) THEN \
            IF (tot_charge >= dft_chg_max) THEN \
               tot_charge = prev_chg2 + 0.7*(dft_chg_max-prev_chg2) \
            ELSE \
               tot_charge = tot_charge + 0.7*(dft_chg_max-tot_charge) \
            END IF \
         ELSE IF ((tot_charge + change_vec) < dft_chg_min) THEN \
            IF (tot_charge <= dft_chg_min) THEN \
               tot_charge = prev_chg2 - 0.7*(prev_chg2-dft_chg_min) \
            ELSE \
               tot_charge = tot_charge - 0.7*(tot_charge-dft_chg_min) \
            END IF \
 \
         ELSE \
            tot_charge = tot_charge + change_vec \
 \
         END IF \
         WRITE(environ_unit,*)"DFT_min ",dft_chg_min \
         WRITE(environ_unit,*)"DFT_max ",dft_chg_max \
         CALL mp_bcast(tot_charge, ionode_id,intra_image_comm) \
         !print *,"DFT_max",dft_chg_max \
         cur_chg = tot_charge \
         prev_step_size = ABS(cur_chg - prev_chg) \
         prev_dchg = cur_dchg \
         WRITE(environ_unit,*)"Convergeable? ",converge \
         CALL mp_bcast(converge,ionode_id, intra_image_comm) \
         CALL mp_bcast(prev_step_size,ionode_id,intra_image_comm) \
         IF (((prev_step_size > chg_thr) .OR. (.NOT. converge)) & \
                        & .AND. (chg_step < nstep-1))  THEN \
             conv_ions = .FALSE. \
             WRITE( STDOUT, 1002)& \
               &chg_step,cur_fermi,ss_chg,prev_step_size,cur_dchg,tot_charge \
             CALL qexsd_set_status(255) \
             CALL punch( "config" ) \
             CALL add_qexsd_step(i_step) \
             nelec = ionic_charge - tot_charge \
             CALL mp_bcast(nelec, ionode_id,intra_image_comm) \
             CALL update_pot() \
             CALL hinit1() \
         ELSE \
             IF (chg_step == nstep -1) THEN \
                WRITE(STDOUT,*)NEW_LINE("a")//"   Exceeded Max number steps!"//& \
                         &NEW_LINE("a")//"   Results probably out of accurate range"//& \
                         &NEW_LINE("a")//"   Smaller chg_thr recommended."//& \
                         &NEW_LINE("a")//"   Writing current step to q-v.dat." \
             END IF \
             WRITE(STDOUT, 1003)chg_step,prev_step_size,ss_chg,cur_dchg,& \
                                &bulk_potential \
             OPEN(21,file = "q-v.dat", status = "unknown") \
             WRITE(37, *)"Potential (V-V_fb)  Surface State Potential (V-V_cut)",& \
                          &"  Electrode Charge (e)",& \
                          &"  Surface States Charge (e)    ",& \
                          &"Electrode Charge per surface area (e/cm^2)     ",& \
                          &"Surface State Charge per surface area (e/cm^2)" \
             sqbohrtosqcm = 2.8002D-17     ! Value of 1 square bohr in cm^2 \
             surf_area = area*sqbohrtosqcm \
             chg_per_area = electrode_charge/surf_area \
             ss_chg_per_area = ss_chg/surf_area \
             ss_potential = -bulk_potential \
             CALL mp_bcast(ss_potential, ionode_id, intra_image_comm) \
             !print *, bulk_potential,ss_potential \
             WRITE(37, 1004)total_potential, ss_potential,&\
                             &electrode_charge, ss_chg,& \
                             &chg_per_area,ss_chg_per_area \
             CLOSE(21) \
         END IF \
    END IF \
 \
    CALL stop_clock( "semiconductor" ) \
  END IF \
 \
 \
 \
1001 FORMAT(5x,//"***************************************************",//& \
                &"     Flatband potential calculated as ",F14.8,// & \
                &"     Now using initial charge of:  ",F14.8,// & \
                 "***************************************************") \
  ! \
1002 FORMAT(5x,//"***************************************************",//& \
                &"     Finished Charge convergence step : ",I3,//& \
                &"     DFT Fermi level calculated as ",F14.8,// & \
                &"     Charge trapped in surface states: ",F14.8," e",//& \
                &"     Charge Accuracy < ",F14.8,// & \
                &"     Difference between bulk and DFT fermi: ",F14.8,//& \
                &"     Now using DFT charge of:  ",F14.8,// & \
                 "***************************************************") \
1003 FORMAT(5x,//"***************************************************",//& \
                &"     Finished charge convergence step : ",I3,//& \
                &"     Convergence of charge with accuracy < ",F14.8," e",// & \
                &"     Charge trapped in surface states: ",F14.8,//& \
                &"     Difference between bulk and DFT fermi: ",F14.8,//& \
                &"     Final Potential: ",F14.8," V", //& \
                &"     Output written to q-v.dat       ",//& \
                 "***************************************************") \
1004 FORMAT(1x,4F14.8,2ES12.5) \
!Environ patch 
' tmp.2 > tmp.1

mv tmp.1 plugin_semiconductor_update.f90

rm tmp.2
