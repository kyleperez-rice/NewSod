!!****if* source/Simulation/SimulationMain/NewSod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockId)
!!                       
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod problem
!!
!!  References:  Hi
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  
!!
!! PARAMETERS
!!
!!  sim
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_minRhoInit     Density floor for initial condition
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)
	! Get the constants
	#include "constants.h"
	#include "Flash.h"
	#include "Eos.h"
	
	! Get the needed unit scope data
	use Simulation_data, ONLY: 	sim_posn, sim_xCos, sim_yCos, sim_zCos,	&
					sim_rhoLeft, sim_pLeft, sim_uLeft,		&
					sim_rhoRight, sim_pRight, sim_uRight,		&
					sim_smallX, sim_gamma, sim_smallP
	use Grid_interfaces, ONLY :	Grid_getBlkIndexLimits, Grid_getCellCoords, &
					Grid_putPointData
	implicit none
	
	
	
	! Define arguments and indicate if they are inputs or outputs
	integer, intent(in) :: blockID
	
	! Declare all local variables
	integer :: i, j, k, n
	integer :: iMax, jMax, kMax
	real :: xx, yy, zz, xxL, xxR
	real :: lPosn0, lPosn
	
	! Arrays to hold coordinate information for a given block
	real, allocatable, dimension(:) :: xCenter, xLeft, xRight, yCoord, zCoord
	
	! Array to get integer indices defining the beginning and the end of a block
	integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
	
	! The number of grid points along each dimension
	integer :: sizeX, sizeY, sizeZ
	
	integer, dimension(MDIM) :: axis
	integer :: dataSize
	logical :: gcell = .true
	
	! These variables store the calculated initial values of physical variables at a certain grid point at some time
	real :: rhoZone, velxZone, velyZone, velzZone, presZone, &
		enerZone, ekinZone
		
	! Get the integer endpoints of the block in all dimensions
	! The Array blkLimits returns the interior end points
	! The array blkLimitsGC returns endpoints, including guardcells
	call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
	
	! Get the size along each dimension for allocation and then allocate
	sizeX = blkLimitsGC(HIGH,IAXIS)
	sizeY = blkLimitsGC(HIGH,JAXIS)
	sizeZ = blkLimitsGC(HIGH,KAXIS)
	allocate(xLeft(sizeX))
	allocate(xRight(sizeX))
	allocate(xCenter(sizeX))
	allocate(yCoord(sizeY))
	allocate(zCoord(sizeZ))
	
	! Setting up initial conditions. In this example, we get the coordinates for the cells in their current block.
	xCoord(:) = 0.0
	yCoord(:) = 0.0
	zCoord(:) = 0.0
	
	call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE, gcell, xLeft, sizeX)
	call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
	call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, gcell, xRight, sizeX)
	call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCoord, sizeY)
	call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCoord, sizeZ)
	
	!---
	! Loop over all the zones in the current block and set the variables.
	!---
	do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
		zz = zCoord(k) ! Coordinates of the cell center in the z-direction
		lPosnO = sim_posn - zz*sim_zCos/sim_xCos	! Where along the x-axis
								! the shock intersects
								! the xz-plane at the current z
		do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
			yy = yCoord(j) ! Center coordinates in the y direction
			lPosn = lPosn0 - yy*sim_yCos/sim_xCos		! The position of the shock in
									! the current yz row
			dataSize = 1	! For Grid put data function, we are initializing
					! a single cell at a time and sending it to Grid
				
			do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
				xx = xCenter(i)	! Center coordinate along x
				xxL = xLeft(i)		! Left coordinate along y
				xxR = xRight(i)	! Right coordinate along z
	! For Sod problem, we have a discontinuity along the shock plane.
	! We do this by initializing the grid points to the left of the shock plane with one value
	! And the points to the right with another value.
	
				if (xxR <= lPosn) then
					rhoZone = sim_rhoLeft
					presZone = sim_pLeft
		
					velxZone = sim_uLeft * sim_xCos
					velyZone = sim_uLeft * sim_yCos
					velzZone = sim_uLeft * sim_zCos
		
		! Initialize the cells which straddle the shock.
		! Treat them as if one half lays to the left and the other half to the right.
				elseif ((xxL < lPosn) .and. (xxR > lPosn)) then
	
					rhoZone = 0.5 * (sim_rhoLeft+sim_rhoRight)
					presZone = 0.5 * (sim_pLeft + sim_pRight)
		
					velxZone = 0.5 *(sim_uLeft+sim_uRight) * sim_xCos
					velyZone = 0.5 *(sim_uLeft+sim_uRight) * sim_yCos
					velzZone = 0.5 *(sim_uLeft+sim_uRight) * sim_zCos
				! Intialize cells to the right of the initial shock.
				else
	
					rhoZone = sim_rhoRight
					presZone = sim_pRight
		
					velxZone = sim_uRight * sim_xCos
					velyZone = sim_uRight * sim_yCos
					velzZone = sim_uRight * sim_zCos
				endif
				! Get the position of the cell in the block
				axis(IAXIS) = i	
				axis(JAXIS) = j
				axis(KAXIS) = k
	
	! Compute the gass energy and set the gamma values needed for the equation of state
	
				ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
	
				enerZone = presZone / (sim_gamma -1.)
				enerZone = enerZone / rhoZone
				enerZone = enerZone + ekinZone
				enerZone = max(enerZone, sim_smallP)
	
	! Store the variables in the current zone via the Grid_putPointData method
				call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
				call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
				call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
				call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
				call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
				call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
				call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
				call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
				enddo
			enddo
		enddo

	deallocate(xLeft)
	deallocate(xRight)
	deallocate(xCenter)
	deallocate(yCoord)
	deallocate(zCoord)


	return
end subroutine Simulation_initBlock(blockID)


