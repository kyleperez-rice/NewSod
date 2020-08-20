!!****if* source/Simulation/SimulationMain/NewSod/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sod
!!  
!! PARAMETERS
!!
!!  sim_rhoLeft	Initial Density to the left
!!  sim_rhoRight	Initial Density to the right
!!  sim_pLeft		Intial pressure...
!!  sim_pRight		Intial pressure...
!!
!!  sim_uLeft		Initial velocity...
!!  sim_uRight		Initial velocity...
!!  sim_xAngle		Angle of the discontinuity with respect to the x axis
!!  sim_yAngle		Angle of the discontinuity with respect to the y axis
!!  sim_posn		The intersection between the shock plane and the x axis
!!
!!  sim_gamma		Ratio of specific heats c_p/c_v
!!  sim_smallP		???
!!  sim_smallX		???
!!
!!  sim_xCos		???
!!  sim_yCos		???	
!!  sim_zCos		???	
!!
!!***

module Simulation_data
#include "Flash.h"
	implicit none
	
	!! Runtime Parameters
	real, save :: sim_rhoLeft, sim_rhoRight, sim_pLeft, sim_pRight
	real, save :: sim_uLeft, sim_uRight, sim_xAngle, sim_yAngle, sim_posn
	real, save :: sim_gamma, sim_smallP, sim_smallX
	
	!! Other unit variables
	real, save :: sim_xCos, sim_yCos, sim_zCos

end module Simulation_data
