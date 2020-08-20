!!****if* source/Simulation/SimulationMain/NewSod/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sod problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!  sim_smallP		???
!!  sim_smallX		???
!!  sim_gamma		Ratio of the specific heats
!!  sim_rhoLeft	Initial density on the left
!!  sim_rhoRight	Initial density on the right
!!  sim_pLeft		Initial pressure...
!!  sim_pRight		Initial pressure...
!!  sim_uLeft		Initial velocity...
!!  sim_uRight		Initial velocity...
!!  sim_xAngle		Angle of the discontinuity with respect to the x axis
!!  sim_yAngle		...with respect to the y axis
!!  sim_posn		Intersection of the shock plane and the x axis
!!
!!***

subroutine Simulation_init()
	
	use Simulation_data
	use RuntimeParameters_interface, ONLY : RuntimeParameters_get
	implicit none
	#include "Flash.h"
	#include "constants.h"
		! get the relevant runtime parameters for the problem
		
		call RuntimeParameters_get('smallp', sim_smallP)
		call RuntimeParameters_get('smallx', sim_smallX)
		call RuntimeParameters_get('gamma', sim_gamma)
		call RuntimeParameters_get('sim_rhoLeft', sim_rhoLeft)
		call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
		call RuntimeParameters_get('sim_pLeft', sim_pLeft)
		call RuntimeParameters_get('sim_pRight', sim_pRight)
		call RuntimeParameters_get('sim_uLeft', sim_uLeft)
		call RuntimeParameters_get('sim_uRight', sim_uRight)
		call RuntimeParameters_get('sim_xangle', sim_xAngle)
		call RuntimeParameters_get('sim_yangle', sim_yAngle)
		call RuntimeParameters_get('sim_posn', sim_posn)
		
		! Do other initializations
		! Convert the shock angle parameters
		sim_xAngle = sim_xAngle * 0.0174532925 ! Converts degrees to radians
		sim_yAngle = sim_yAngle * 0.0174532925 ! Converts degrees to radians
		
		sim_xCos = cos(sim_xAngle)
		
		if (NDIM == 1) then
			sim_xCos = 1.
			sim_yCos = 0.
			sim_zCos = 0.
			
		elseif (NDIM == 2) then
			sim_yCos = sqrt(1. -sim_xCos*sim_xCos)
			sim_zCos = 0.
			
		elseif (NDIM == 3) then
			sim_yCos = cos(sim_yAngle)
			sim_zCos = sqrt( max(0., 1. -sim_xCos*sim_xCos - sim_yCos*sim_yCos) )
		endif
end subroutine Simulation_init
