!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_semiconductor_update(i_step)
  !----------------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout, ionode
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir
  USE io_global,      ONLY: ionode, ionode_id
  USE mp,             ONLY: mp_bcast, mp_barrier, mp_sum
  USE mp_world,       ONLY: world_comm
  USE mp_images,      ONLY: intra_image_comm
  USE mp_bands,       ONLY: intra_bgrp_comm
  USE klist,            ONLY : tot_charge, nelec
  USE cell_base,        ONLY : omega
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : rho
  USE control_flags,    ONLY : conv_ions, nstep
  USE ener,             ONLY : ef
  USE constants,        ONLY : rytoev
  USE fft_base,         ONLY : dfftp
  USE ions_base,        ONLY : nat, ityp, zv
  USE extrapolation,    ONLY : update_pot
  !
  USE plugin_flags
  USE qexsd_module,     ONLY:   qexsd_set_status


! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
  


  !
  IMPLICIT NONE
  !

  SAVE 




  INTEGER,INTENT(INOUT)     ::   i_step

!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***

! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!


END SUBROUTINE plugin_semiconductor_update

