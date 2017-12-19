!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_initialization()
  !----------------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout, ionode
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir
  USE klist,            ONLY : tot_charge
  USE control_flags,    ONLY : nstep
  !
  USE plugin_flags
  !
  
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!

  IMPLICIT NONE
  !
  !
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***

END SUBROUTINE plugin_initialization