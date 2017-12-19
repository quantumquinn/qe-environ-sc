! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by O. Andreussi and N. Marzari



!
! this module is designed to find the electrochemical response of semiconductor
! electrodes in contact with water.
!
! original by Q. Campbell qjc5019@psu.edu
!



!
!--------------------------------------------------------------------
MODULE semiconductor 
!--------------------------------------------------------------------

  USE kinds,          ONLY: DP
  USE constants,      ONLY: pi, e2
  USE io_global,      ONLY: stdout
  USE mp,             ONLY: mp_sum
  USE mp_bands,       ONLY: intra_bgrp_comm
  USE environ_cell,   ONLY: domega, omega, at, ibrav
  USE environ_ions,   ONLY: avg_pos, rhoions
  USE environ_base,   ONLY: verbose, environ_unit,                   &
                            env_external_charges, extcharge_dim,     &
                            extcharge_axis, extcharge_pos,           &
                            extcharge_spread, extcharge_charge,      &
                            rhoexternal, semiconductor_region,       &
                            slab_direction, helmholtz_distance,      &
                            electrode_charge, sc_dielectric,         &
                            env_dielectric_regions, epsregion_eps,   &
                            epsregion_dim, epsregion_axis,           &
                            epsregion_pos, epsregion_width,          &
                            epsregion_spread, sc_cutoff,             &
                            dopant_concentration, flatband_fermi,    &
                            macro_avg_window                           
  !


  IMPLICIT NONE

  SAVE 

 
  INTEGER :: slab_dir, discont_loc
  !
  REAL(DP) :: max_atomic_pos, min_atomic_pos
  REAL(DP), ALLOCATABLE :: macro_pot_0(:)
  REAL(DP), ALLOCATABLE :: current_v(:)
  REAL(DP), ALLOCATABLE :: z_0(:)
  REAL(DP)                :: bulk_potential
  REAL(DP)                :: area
  !
  LOGICAL :: insert_helmholtz  

  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE check_slab(tau,alat,nat, at)
!--------------------------------------------------------------------
  USE kinds,          ONLY: DP

  IMPLICIT NONE

  REAL(DP) :: slab_dir_distance 
  INTEGER, INTENT( IN )     :: nat
  REAL(DP), INTENT( IN ) :: tau(3,nat)
  REAL ( DP ), INTENT( IN ) :: alat
  REAL( DP ), INTENT( IN ) :: at(3,3)

  !Assigning numerical slab direction
  SELECT CASE( TRIM(slab_direction))

  CASE('x')
    !
    slab_dir =1
    area=at(2,2)*at(3,3)*alat**2 
  CASE('y')
    slab_dir = 2
    area=at(1,1)*at(3,3)*alat**2
  CASE('z')
    slab_dir = 3
    area=at(1,1)*at(2,2)*alat**2
  CASE DEFAULT
    CALL errore('semiconductor','slab_direction = ' // &
                & trim(slab_direction) // ' not recognized',1)


  END SELECT


  !
  !Finds the maximum atomic value in the slab direction
  !If this value is not > 3 angstrom from the cell edge,
  !it will throw an error
  !
  max_atomic_pos = MAXVAL(tau(slab_dir, :)*alat)
  min_atomic_pos = MINVAL(tau(slab_dir, :)*alat)
  IF (ibrav == 0) THEN
     slab_dir_distance = alat * at(slab_dir,slab_dir)
     IF ((slab_dir_distance - max_atomic_pos) < helmholtz_distance .OR. &
         &(min_atomic_pos < helmholtz_distance )) THEN
       CALL errore('semiconductor','helmholtz distance would put external charge'// &
                    & NEW_LINE('C') //'     outside of cell. Recheck dimensions!',1)
     END IF 
     IF (min_atomic_pos < 5.669) THEN !Using a 3 angstrom distance as min(?)
       CALL errore('semiconductor','Not 3 angstrom of vaccuum on sc side'// &
                     & NEW_LINE('C') //'     Recheck dimensions!',2)
     END IF


  ELSE
    CALL errore('semiconductor','Only ibrav=0 implemented',1)
  END IF


  !WRITE(environ_unit, 200)

  RETURN

200   FORMAT(1X,'Input is a slab...')
!--------------------------------------------------------------------
  END SUBROUTINE check_slab
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE generate_sc_extcharge(tot_charge)
!--------------------------------------------------------------------
  USE kinds,          ONLY: DP

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: tot_charge  ! Charge on DFT slab

  REAL(DP) :: slab_z

  !
  !  Creates planar Helmholtz counter charge on the soln side
  !  and a matching one on the semiconductor side to ensure charge
  !  neutrality
  !


  IF (insert_helmholtz) THEN
     slab_z = max_atomic_pos + helmholtz_distance
     WRITE( environ_unit, * )&
                &"slab distance: ",slab_z
     !extcharge_pos(1,env_external_charges) = 0.D0
     !extcharge_pos(2,env_external_charges) = 0.D0
     !extcharge_pos(3,env_external_charges) = 0.D0
     extcharge_pos(slab_dir,env_external_charges-1) = slab_z
     !WRITE( environ_unit, *)'external charge pos: ',extcharge_pos
     extcharge_charge(env_external_charges-1) = electrode_charge
     extcharge_spread(env_external_charges-1) = 0.5
     extcharge_dim(env_external_charges-1) = 2 ! Assuming planar
     extcharge_axis(env_external_charges-1) = slab_dir

     slab_z = min_atomic_pos - helmholtz_distance
     extcharge_pos(slab_dir,env_external_charges) = slab_z
     !WRITE( environ_unit, *)'external charge pos: ',extcharge_pos
     extcharge_charge(env_external_charges) = tot_charge-electrode_charge
     extcharge_spread(env_external_charges) = 0.5
     extcharge_dim(env_external_charges) = 2 ! Assuming planar
     extcharge_axis(env_external_charges) = slab_dir

  END IF 

  RETURN


!--------------------------------------------------------------------
  END SUBROUTINE generate_sc_extcharge
!--------------------------------------------------------------------


!--------------------------------------------------------------------
  SUBROUTINE generate_sc_epsregion()
!--------------------------------------------------------------------
  USE kinds,          ONLY: DP

  IMPLICIT NONE

  REAL(DP) :: half_sc, half_slab
  !
  !Creates two dielectric regions for the slab:
  !  one of the sc dielectric on the sc vacuum side
  !  and one of dielectric 1 within the slab to help
  !  convergence
  !
  half_sc = min_atomic_pos/2.D0
  half_slab = (max_atomic_pos - min_atomic_pos)/2.D0

  WRITE(environ_unit, *)&
               &"Dielectric regions: ",env_dielectric_regions

  !Making first dielectric region on sc side


  epsregion_eps(1,env_dielectric_regions-1) = sc_dielectric
  epsregion_eps(2,env_dielectric_regions-1) = sc_dielectric
  epsregion_pos(slab_dir,env_dielectric_regions-1) = half_sc
  epsregion_width(env_dielectric_regions-1) = half_sc - 0.125_DP
  epsregion_spread(env_dielectric_regions-1) = 0.125 !Assuming small spread
  epsregion_dim(env_dielectric_regions-1) = 2 !Making planar section
  epsregion_axis(env_dielectric_regions-1) = slab_dir
  
  !Making second dielectric region within slab

  epsregion_eps(1,env_dielectric_regions) = 1.0_DP
  epsregion_eps(2,env_dielectric_regions) = 1.0_DP
  epsregion_pos(slab_dir,env_dielectric_regions) = min_atomic_pos + half_slab
  epsregion_width(env_dielectric_regions) = half_slab -2.0_DP
  epsregion_spread(env_dielectric_regions) = 0.25 !Assuming small spread
  epsregion_dim(env_dielectric_regions) = 2 !Making planar section
  epsregion_axis(env_dielectric_regions) = slab_dir

!--------------------------------------------------------------------
  END SUBROUTINE generate_sc_epsregion
!--------------------------------------------------------------------

! This saves the vfull_free at each step to the variable current_v
! which can be used anywhere in this module (and will be)
! seems inefficient but here we go

!--------------------------------------------------------------------
  SUBROUTINE save_current_v(nnr,v)
!--------------------------------------------------------------------


  USE environ_base,   ONLY: environ_unit
  USE io_global,      ONLY: ionode, ionode_id
  USE mp,             ONLY: mp_bcast
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_images,      ONLY: intra_image_comm
  USE fft_base,       ONLY: dfftp
  USE scatter_mod,    ONLY: gather_grid

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nnr
  REAL(DP), INTENT(IN)   :: v(nnr)
  IF ( ALLOCATED( current_v)) DEALLOCATE(current_v)
  ALLOCATE(current_v(nnr))
  current_v = 0.0
  current_v = v
  WRITE(environ_unit,*)"Saved vfull_free"

  RETURN
!--------------------------------------------------------------------
  END SUBROUTINE save_current_v
!--------------------------------------------------------------------




! This routing saves the flatband pot to the variable vltot_zero
! which can be used anywhere in this module (and will be)

!--------------------------------------------------------------------
  SUBROUTINE save_flatband_pot(nnr)
!--------------------------------------------------------------------
 ! Used to analyze the potential at the end of flatband calculation and to save
 ! the final results in semiconductor variables 


  USE environ_base,   ONLY: environ_unit
  USE kinds,          ONLY: DP
  USE io_global,      ONLY: ionode, ionode_id
  USE mp,             ONLY: mp_bcast, mp_barrier, mp_gather
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_world,       ONLY: world_comm 
  USE mp_images,      ONLY: intra_image_comm
  USE fft_base,       ONLY: dfftp
  USE scatter_mod,    ONLY: gather_grid

  IMPLICIT NONE


  INTEGER, INTENT(IN)   :: nnr
  REAL(DP)  :: vltot_zero(nnr)
  REAL(DP), ALLOCATABLE :: z(:)
  REAL(DP), ALLOCATABLE :: planar_v(:)
  REAL(DP), ALLOCATABLE :: macro_v(:)
  INTEGER               :: nnrz
  REAL(DP)              :: avg_window


! Should have user option to set avg window eventually
  avg_window = macro_avg_window
  vltot_zero = current_v 
  DEALLOCATE(current_v)
  nnrz = 0
  CALL find_planar_avg(nnr,vltot_zero,z,planar_v,nnrz,0)
  ALLOCATE(macro_v(nnrz))
  CALL fix_discont(nnrz,z,planar_v,0)
  CALL align_potential(nnrz,z,planar_v,flatband_fermi)
  CALL find_macro_avg(nnrz,z,planar_v,macro_v,avg_window,0)
  ALLOCATE(macro_pot_0(nnrz))
  ALLOCATE(z_0(nnrz))
  z_0 = z
  macro_pot_0 = macro_v

  CALL mp_bcast(flatband_fermi,  ionode_id, intra_image_comm)


  WRITE(environ_unit,*)"Saved flatband potential"


  !CALL mp_barrier(world_comm)

  RETURN
!--------------------------------------------------------------------
  END SUBROUTINE save_flatband_pot
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE read_flatband_pot()
!--------------------------------------------------------------------
 ! Used to analyze the potential at the end of flatband calculation and to save
 ! the final results in semiconductor variables 

  USE environ_base,   ONLY: environ_unit
  USE kinds,          ONLY: DP
  USE io_global,      ONLY: ionode, ionode_id
  USE mp,             ONLY: mp_bcast, mp_barrier, mp_gather
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_world,       ONLY: world_comm
  USE mp_images,      ONLY: intra_image_comm
  USE fft_base,       ONLY: dfftp
  USE scatter_mod,    ONLY: gather_grid
  USE environ_base,   ONLY: flatband_dir

  IMPLICIT NONE


  REAL(DP), ALLOCATABLE :: z(:)
  REAL(DP), ALLOCATABLE :: planar_v(:)
  REAL(DP), ALLOCATABLE :: macro_v(:)
  INTEGER               :: nnrz
  REAL(DP)              :: avg_window
  INTEGER               :: old_planar_unit
  INTEGER               :: ios, i
  LOGICAL               :: ext
  REAL(DP)              :: dummy_fermi

! Should have user option to set avg window eventually
  avg_window = macro_avg_window
  
  old_planar_unit = 33
  nnrz = 0
  ios = 0
  INQUIRE(file = TRIM(flatband_dir)//"planar_pot_000.dat", exist = ext)

  IF( .NOT. ext ) CALL errore('read_flatband', &
                      & "can't find planar_pot_000.dat", 1)

  OPEN(old_planar_unit, file=TRIM(flatband_dir)//"planar_pot_000.dat")
  DO 
    READ(old_planar_unit,*, iostat = ios)
    IF (ios /=0) EXIT
    nnrz = nnrz + 1
  END DO

  REWIND(old_planar_unit)

  ALLOCATE(z(nnrz))
  ALLOCATE(planar_v(nnrz))

  DO i = 1, nnrz
    READ(old_planar_unit,FMT=1002)z(i),planar_v(i)
  END DO
  CLOSE(old_planar_unit)
  ALLOCATE(macro_v(nnrz))
  dummy_fermi = 0.0
  CALL fix_discont(nnrz,z,planar_v,0)
  CALL align_potential(nnrz,z,planar_v,dummy_fermi)
  CALL find_macro_avg(nnrz,z,planar_v,macro_v,avg_window,0)
  ALLOCATE(macro_pot_0(nnrz))
  ALLOCATE(z_0(nnrz))
  z_0 = z
  macro_pot_0 = macro_v


  !WRITE(environ_unit,*)"Read flatband potential"


  !CALL mp_barrier(world_comm)

  RETURN

1002  FORMAT(2(2X,2f14.10))
!--------------------------------------------------------------------
  END SUBROUTINE read_flatband_pot
!--------------------------------------------------------------------


!--------------------------------------------------------------------
  SUBROUTINE save_current_pot(nnr,fermi_q,dq,ss_chg,v_cut,istep)
!--------------------------------------------------------------------

! Used to save and find the "derivative" of charge for scf runs once charge
! starts to be varied works the same as save_flatband_pot




  USE environ_base,   ONLY: environ_unit
  USE kinds,          ONLY: DP
  USE io_global,      ONLY: ionode, ionode_id
  USE mp,             ONLY: mp_bcast, mp_barrier, mp_gather
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_world,       ONLY: world_comm
  USE mp_images,      ONLY: intra_image_comm
  USE fft_base,       ONLY: dfftp
  USE scatter_mod,    ONLY: gather_grid

  IMPLICIT NONE

  INTEGER,  INTENT(IN)    :: nnr, istep
  REAL(DP), INTENT(INOUT)  :: fermi_q 
  REAL(DP), INTENT(INOUT)  :: dq
  REAL(DP), INTENT(OUT)   :: v_cut
  REAL(DP), INTENT(OUT) :: ss_chg
  REAL(DP), ALLOCATABLE :: vltot_zero(:)
  REAL(DP), ALLOCATABLE :: z(:)
  REAL(DP), ALLOCATABLE :: planar_v(:)
  REAL(DP), ALLOCATABLE :: macro_v(:)
  INTEGER               :: nnrz
  REAL(DP)              :: avg_window


! Should have user option to set avg window eventually

  avg_window = macro_avg_window
  ALLOCATE(vltot_zero(nnr))
  vltot_zero = current_v
  DEALLOCATE(current_v)
  nnrz = 0
  CALL find_planar_avg(nnr,vltot_zero,z,planar_v,nnrz,istep)
  ALLOCATE(macro_v(nnrz))
  CALL fix_discont(nnrz,z,planar_v,istep)
  CALL align_potential(nnrz,z,planar_v,fermi_q)
  CALL find_macro_avg(nnrz,z,planar_v,macro_v,avg_window,istep)
  CALL d_chg(nnrz, z,macro_v, fermi_q, dq, ss_chg,v_cut,istep)

  WRITE(environ_unit,*)"Potential Analyzed"
  CALL mp_bcast(fermi_q,  ionode_id, intra_image_comm)
  CALL mp_bcast(dq,       ionode_id, intra_image_comm)
  CALL mp_bcast(ss_chg,   ionode_id, intra_image_comm)

  CALL mp_barrier(world_comm)

  RETURN
!--------------------------------------------------------------------
  END SUBROUTINE save_current_pot
!--------------------------------------------------------------------


! This method will find the planar average of a given
! qe fft grid based file

!--------------------------------------------------------------------
  SUBROUTINE find_planar_avg(nnr,pot_file,z,pot_avg,nnrz,istep)
!--------------------------------------------------------------------
  USE environ_base,   ONLY: environ_unit
  USE fft_base,       ONLY: dfftp
  USE cell_base,      ONLY: at, alat
  USE io_global,      ONLY: ionode
  USE kinds,          ONLY: DP
! Compatible with QE-5.3 svn
  USE scatter_mod,    ONLY : gather_grid
  USE mp,             ONLY : mp_sum
  USE mp_bands,       ONLY : intra_bgrp_comm

  IMPLICIT NONE 

  INTEGER, INTENT(IN)    :: nnr
  INTEGER, INTENT(OUT)   :: nnrz
  INTEGER, INTENT(IN)    :: istep
  REAL(DP), INTENT(IN)   :: pot_file(nnr)
  ! copying variables from write cube method
  ! Some may not be necessary
  INTEGER                  :: nr1, nr2, nr3
  INTEGER                  :: nr1x, nr2x, nr3x
  INTEGER                  :: ir, i, j, k, planar_unit
      !
  REAL( DP ), INTENT(OUT), ALLOCATABLE  :: pot_avg( : )
  REAL( DP ), INTENT(OUT), ALLOCATABLE  :: z( : )
  REAL( DP ), ALLOCATABLE  :: flocal( : )
  REAL( DP ), DIMENSION(3) :: origin
  CHARACTER (LEN=8)        :: fmt
  CHARACTER (LEN=3)        :: x1
      !
      !
  nr1x = dfftp%nr1x
  nr2x = dfftp%nr2x
  nr3x = dfftp%nr3x


  nr1 = dfftp%nr1
  nr2 = dfftp%nr2
  nr3 = dfftp%nr3
      !
  ALLOCATE( flocal( nr1x*nr2x*nr3x ) )
  flocal = 0.D0
 !Compatible with QE-svn 5.3
! I'm just copying these lines from the write_cube method
! I suspect they are mainly distributing the points across
! processors. But really I don't know nothing
  CALL gather_grid( dfftp,pot_file, flocal )
  CALL mp_sum( flocal, intra_bgrp_comm )


  !flocal = pot_file
  origin = 0.D0

  planar_unit = 700


  WRITE(environ_unit, 400)

! Finds the macroscopic avg based on the direction
! Asssumes that slab_dir has already been checked and can only
! be 1,2,3

! This can only work for one processor
  IF( ionode ) THEN
     fmt = '(I3.3)'
     write(x1,fmt)istep      
     IF (( verbose .GE. 3) .OR. (istep == 0))&
      &OPEN( planar_unit, file = 'planar_pot_'//trim(x1)//'.dat',&
                & status = 'unknown')
     IF ( slab_dir ==1) THEN
        ALLOCATE( pot_avg(nr1))
        ALLOCATE( z(nr1))
        nnrz = nr1
        DO i = 1, nr1
           pot_avg(i) = 0.D0
           DO j = 1 , nr2
               DO k = 1, nr3
                  ir = i + (j-1) *nr1 + (k - 1)*nr1*nr2
                  pot_avg(i) = pot_avg(i) + dble(flocal(ir))
               ENDDO        
           ENDDO
           pot_avg(i) = pot_avg(i)/(dble(nr2*nr3))
           z(i) = origin(1) + i*at(1,1)/dble(nr1)*alat
           IF (( verbose .GE. 3 ) .OR. (istep == 0))&
                         & WRITE(planar_unit,401),z(i),pot_avg(i) 
        ENDDO
     ELSEIF (slab_dir == 2) THEN
        ALLOCATE( pot_avg(nr2))
        ALLOCATE( z(nr2))
        nnrz = nr2
        DO j =1, nr2
           pot_avg(j) = 0.D0
           DO i =1, nr1
              DO k =1, nr3
                 ir = i + (j -1)*nr1 + (k-1)*nr1*nr2
                 pot_avg(j) = pot_avg(j) + dble(flocal(ir))
              ENDDO
           ENDDO
           pot_avg(j) = pot_avg(j)/(dble(nr1*nr3))
           z(j) = origin(2) + j*at(2,2)/dble(nr2)*alat
           IF (( verbose .GE. 3) .OR. (istep == 0))&
                   &WRITE(planar_unit,401),z(j),pot_avg(j)
        ENDDO
     ELSE
       !WRITE(planar_unit,*)"Checking nr1,nr2,nr3",nr1,nr2,nr3
       !WRITE(planar_unit,*)"Checking nr1x, nr2x",nr1x,nr2x,nr3x
       ALLOCATE( pot_avg(nr3))
       ALLOCATE( z(nr3))
       nnrz = nr3
       DO k =1, nr3
           pot_avg(k) = 0.D0
           DO j =1, nr2
              DO i =1, nr1
                 ir = i + (j -1)*nr1 + (k-1)*nr1*nr2
                 pot_avg(k) = pot_avg(k) + dble(flocal(ir))
                 !WRITE(environ_unit,*)"Pot_avg = ",pot_avg(k)
              ENDDO
           ENDDO
           pot_avg(k) = pot_avg(k)/( dble( nr1*nr2 ) )
           !WRITE(planar_unit,*)'nr1*nr2' ,( dble( nr1*nr2 ) )
           z(k) = origin(3) + k*at(3,3)/dble(nr3)*alat
           IF (( verbose .GE. 3) .OR. (istep == 0))&
                        &WRITE(planar_unit,401),z(k),pot_avg(k)
        ENDDO
     END IF
     CLOSE(planar_unit)
  END IF

  DEALLOCATE( flocal )

  RETURN
400   FORMAT(1X,'Finding Planar average')
401   FORMAT(2(2X,2f14.10))
!--------------------------------------------------------------------
  END SUBROUTINE find_planar_avg
!--------------------------------------------------------------------

!
! Hopefully by having this all in the same program it will avoid
! the issues of mismatched macro avgs that we saw earlier
!

!--------------------------------------------------------------------
  SUBROUTINE find_macro_avg(nnrz, z,pot_avg,macro_pot,avg_window, istep)
!--------------------------------------------------------------------

! Used to find the macroscopic avg of a planar potentail
! Assumes z and pot already separated and have size nnrz (multiprocess might
! still be an issue on this, unsure)

  USE environ_base,   ONLY: environ_unit
  USE io_global,      ONLY: ionode
  USE kinds,          ONLY: DP

  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: nnrz
  INTEGER, INTENT(IN)     :: istep
  REAL(DP), INTENT(IN)    :: z(nnrz)
  REAL(DP), INTENT(IN)    :: pot_avg(nnrz)
  REAL(DP), INTENT(OUT)   :: macro_pot(nnrz)
  REAL(DP), INTENT(IN)    :: avg_window
  INTEGER                 :: array_window, i, lower, upper
  INTEGER                 :: macro_unit
  CHARACTER (LEN=8)        :: fmt
  CHARACTER (LEN=3)        :: x1

  macro_pot = 0.D0

  macro_unit = 19

  IF( ionode ) THEN
     fmt = '(I3.3)'
     write(x1,fmt)istep
     IF (( verbose .GE. 3) .OR. (istep == 0))&
        &OPEN(macro_unit, file = 'macro_pot_'//trim(x1)//'.dat', status = 'unknown')
     CALL find_z_loc(nnrz,z,avg_window,array_window)
     !WRITE(environ_unit,500)

     IF (MOD(array_window,2) == 1) THEN
        array_window = array_window + 1
     END IF

     DO i=1,nnrz
        lower = i - array_window/2
        upper = i + array_window/2

        IF ( lower < 1) lower =1 
        IF ( upper > nnrz) upper = nnrz
        !WRITE(macro_unit,*)"upper-lower",dble(upper-lower+1)
        !WRITE(macro_unit,*)"pot_avg",pot_avg(lower:upper)
        !WRITE(macro_unit,*)"sum",SUM(pot_avg(lower:upper))
        macro_pot(i) = SUM(pot_avg(lower:upper))/dble(upper-lower+1)
        IF (( verbose .GE. 3) .OR. (istep == 0))& 
                      & WRITE(macro_unit,501)z(i),macro_pot(i)
     END DO
  END IF
  RETURN
500   FORMAT(1X,'Finding Macroscopic Planar average')
501   FORMAT(2(2X,2f14.10))

!--------------------------------------------------------------------
  END SUBROUTINE find_macro_avg
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE find_z_loc(nnrz, z,target_z,i_loc)
!--------------------------------------------------------------------

  !Used to find the array location that matches up to a specific z within the
  !array looking at z. Mainly used to locate derivatives and cutoffs within
  !space

  USE kinds,          ONLY: DP

  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: nnrz
  REAL(DP), INTENT(IN)    :: z(nnrz)
  REAL(DP), INTENT(IN)    :: target_z
  INTEGER, INTENT(OUT)    :: i_loc
  INTEGER                 :: I

  DO i=1,nnrz
     IF (z(i) > target_z) THEN
        i_loc = i
        EXIT
     END IF

  END DO


  RETURN

!--------------------------------------------------------------------
  END SUBROUTINE find_z_loc
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE align_potential(nnrz,z,v,fermi)
!--------------------------------------------------------------------
  ! This aligns the planar potential such that the right edge will be at zero
  USE io_global,      ONLY: ionode

  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: nnrz
  REAL(DP),  INTENT(INOUT) :: v(nnrz)
  REAL(DP),  INTENT(IN) :: z(nnrz)  
  REAL(DP)                :: shift
  REAL(DP), INTENT(INOUT) :: fermi

  ! Currently will shift in to a spot 8 bohr in
  ! probably want it to be within averaging window/2?
  ! might want to think of a more systematic way to do this

  shift = 0.0D0
  IF( ionode) THEN


     shift = v(nnrz)
     v = v - shift

     fermi = fermi - shift

     !WRITE(environ_unit,*)"using shift of ",shift," V"
     !WRITE(environ_unit,*)"adjusted Fermi: ",fermi
  END IF

  RETURN 


!--------------------------------------------------------------------
  END SUBROUTINE align_potential
!--------------------------------------------------------------------
!--------------------------------------------------------------------
  SUBROUTINE fix_discont(nnrz,z,v,istep) 
!--------------------------------------------------------------------
  ! This fixes the apparent discontinuity generated in the electrostatic
  ! potential and ensures that the right edge will be at zero
  USE io_global,      ONLY: ionode
  USE environ_base,   ONLY: environ_unit

  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: nnrz, istep
  REAL(DP),  INTENT(INOUT) :: v(nnrz)
  REAL(DP),  INTENT(IN) :: z(nnrz)
  REAL(DP)            :: v_prev_v, prev_v, cur_v
  REAL(DP)            :: prev_dv, cur_dv
  INTEGER             :: i

  v_prev_v = 0.0
  prev_v = 0.0

  IF( ionode) THEN
     IF (istep == 0) THEN
        discont_loc = 1
        DO i =2, nnrz
           prev_v = v(i-1)
           cur_v = v(i)
           IF (((prev_v < 1.0d-6) .AND. (cur_v > 1.0d-6)) .OR. &
               &((prev_v > 1.0d-6) .AND. (cur_v < 1.0d-6))) THEN
              cur_dv = cur_v - prev_v
              prev_dv = prev_v- v_prev_v
              IF (ABS(cur_dv) > 1.5*ABS(prev_dv)) THEN
                 discont_loc = i
                 EXIT
              END IF 
           END IF
           v_prev_v = prev_v
          
        END DO
     END IF

     !WRITE(environ_unit,*)"Discont_loc : ",discont_loc

     IF ( discont_loc < nnrz/2) THEN
        cur_v = v(discont_loc)
        DO i = 1, discont_loc
           v(i) = cur_v
        END DO
     ELSE
        cur_v = v(discont_loc-1)
        DO i = discont_loc, nnrz
           v(i) = cur_v 
        END DO
     END IF 

  END IF

  RETURN


!--------------------------------------------------------------------
  END SUBROUTINE fix_discont
!--------------------------------------------------------------------



!--------------------------------------------------------------------
  SUBROUTINE d_chg(nnrz, z, v, fermi_q, d_out,ss_chg,v_cut,istep)
!--------------------------------------------------------------------
! Finds how well the fermi levels are aligned. I use this as the derivative of
! charge for the purposes of minimization. Want this to be zero

  USE constants,      ONLY : RYTOEV
  USE io_global,      ONLY: ionode

  IMPLICIT NONE

  INTEGER, INTENT(IN)     :: nnrz
  REAL(DP), INTENT(IN)    :: z(nnrz)
  REAL(DP), INTENT(IN)    :: v(nnrz)
  REAL(DP), INTENT(IN)    :: fermi_q
  REAL(DP), INTENT(OUT)   :: d_out
  REAL(DP), INTENT(OUT) :: ss_chg
  INTEGER, INTENT(IN)     :: istep

  INTEGER                 :: sc_cut_i, i
  REAL(DP), INTENT(OUT)   :: v_cut
  REAL(DP)                :: dv_cut
  REAL(DP)                :: vacuum_perm
  REAL(DP)                :: bohrtoang
  REAL(DP)                :: macro_sub(nnrz)
  REAL(DP)                :: z_ang(nnrz)
  INTEGER                 :: macro_sub_unit
  CHARACTER (LEN=8)        :: fmt
  CHARACTER (LEN=3)        :: x1


  macro_sub_unit = 35

  vacuum_perm = 5.526D21     !This value comes from 
                             ! epsilon_o *(V/angstrom)^2 /(e_o * 1 cm^-3) in V

  bohrtoang = 0.529177

  IF (ionode) THEN


     CALL find_z_loc(nnrz, z, sc_cutoff, sc_cut_i)

     z_ang = z*bohrtoang


     WRITE(environ_unit, *) "sc_cut_loc ", sc_cut_i

     macro_sub = v - macro_pot_0

     v_cut = macro_sub(sc_cut_i) * RYTOEV

     dv_cut = RYTOEV* (macro_sub(sc_cut_i)-macro_sub(sc_cut_i-1))/ &
           & ((z_ang(sc_cut_i) - z_ang(sc_cut_i-1)))  

     WRITE(environ_unit, *) "dv_cut: ",dv_cut

     IF( dv_cut <= 0 ) THEN
       bulk_potential = v_cut + vacuum_perm/2.0 * (dv_cut)**2 /&
                    & dopant_concentration
     ELSE
       bulk_potential = v_cut - vacuum_perm/2.0 * (dv_cut)**2 /&
                    & dopant_concentration
     END IF

     WRITE(environ_unit, *) "bulk_potential: ",bulk_potential


     d_out = -fermi_q +  bulk_potential + flatband_fermi 


     ss_chg = electrode_charge - 0.00154*dv_cut*area !0.00154 comes from 
                                                  !epsilon_0*(V/angstrom)*(1
                                                  !bohr)^2 in units of e

     fmt = '(I3.3)'
     write(x1,fmt)istep

     IF (( verbose .GE. 3) .OR. (istep == 0)) THEN

        OPEN(macro_sub_unit, file = 'macro_subtracted_pot_'//trim(x1)//'.dat', &
                                &status = 'unknown')
        DO i=1,nnrz
           WRITE(macro_sub_unit,701)z(i),macro_sub(i)
        END DO
        CLOSE(macro_sub_unit)
     END IF 
  END IF 

  RETURN
701   FORMAT(2(2X,2f14.10))
!--------------------------------------------------------------------
  END SUBROUTINE d_chg
!--------------------------------------------------------------------




!--------------------------------------------------------------------
END MODULE semiconductor
!--------------------------------------------------------------------

