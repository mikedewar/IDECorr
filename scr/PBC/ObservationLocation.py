from __future__ import division
import pylab as pb
import numpy as pb

def observation_locns(spacestep,estimation_field_width,Delta_s):
	'''Define the center of sensors along x and y
		
		Atguments:
		----------
			spacestep: the spatial step in the simulated field
			observedwidthfield: the width of the observed field
			Delta_s: distance between sensors in mm

		Returns
		-------
			observation_locns_mm:
				the observation location along x or y directions in mm'''
		

	steps_in_field = (2*estimation_field_width)/spacestep + 1;
	inv_spacestep = 1./spacestep;						
	Nspacestep_in_observed_field = inv_spacestep*estimation_field_width+1	

	observation_offest = estimation_field_width/2;     # mm
	observation_offset_units = observation_offest / spacestep -1;
	field_space=pb.arange(-estimation_field_width,estimation_field_width+spacestep,spacestep)
	spatial_location_num=(len(field_space))**2
	Delta_s_units = Delta_s/spacestep	
	nonsymmetric_obs_location_units = pb.arange(1,Nspacestep_in_observed_field,Delta_s_units)
	offset = ((Nspacestep_in_observed_field - nonsymmetric_obs_location_units[-1])/2.)
	symmetricobslocation_units = nonsymmetric_obs_location_units + offset + observation_offset_units

	observation_locs_mm = symmetricobslocation_units*spacestep - estimation_field_width
	return observation_locs_mm

