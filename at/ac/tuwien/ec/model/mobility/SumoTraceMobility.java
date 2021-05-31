package at.ac.tuwien.ec.model.mobility;

import java.io.Serializable;
import java.util.ArrayList;

import at.ac.tuwien.ec.model.Coordinates;

public class SumoTraceMobility implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -546869441881275833L;
	private ArrayList<Coordinates> trace;
	
	public SumoTraceMobility(ArrayList<Coordinates> trace)
	{
		this.trace = trace;
	} 

	public Coordinates getCoordinatesForTimestep(double timestep) {
		int index = (int) Math.ceil(timestep);
		index = index % trace.size();
		return trace.get(index);
	}

}
