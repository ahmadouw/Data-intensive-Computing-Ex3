package at.ac.tuwien.ec.provisioning.edge.mo;

import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.uma.jmetal.operator.mutation.MutationOperator;

public class EdgePlanningMutationOperator implements MutationOperator<EdgePlanningSolution> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7503935801560571134L;
	private double mutationProbability;
	
	public EdgePlanningMutationOperator(double mutationProbability)
	{
		this.mutationProbability = mutationProbability;
	}
	
	@Override
	public EdgePlanningSolution execute(EdgePlanningSolution arg0) {
		UniformRealDistribution mutationRandomGenerator = new UniformRealDistribution();
		UniformIntegerDistribution indexGenerator = 
				new UniformIntegerDistribution(0, arg0.getNumberOfVariables());
		if (mutationRandomGenerator.sample() < mutationProbability)
		{
			int i = indexGenerator.sample();
			Boolean b = arg0.getVariableValue(i);
			arg0.setVariableValue(i, !b);
		}	
		return arg0;
	}

	@Override
	public double getMutationProbability() {
		// TODO Auto-generated method stub
		return 0;
	}

}
