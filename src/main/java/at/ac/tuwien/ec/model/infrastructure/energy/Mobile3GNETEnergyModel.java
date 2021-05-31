package at.ac.tuwien.ec.model.infrastructure.energy;

import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;

import java.io.Serializable;

import at.ac.tuwien.ec.blockchain.Transaction;
import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.SoftwareComponent;


public class Mobile3GNETEnergyModel implements NETEnergyModel,Serializable {

	private final double alpha = 0.025e-3;
	private final double C = 3.5e-3;
	
	
	/* We assume that inData and outData are in bytes
	 * (non-Javadoc)
	 * @see at.ac.tuwien.ec.infrastructuremodel.energy.NETEnergyModel#computeNETEnergy(di.unipi.socc.fogtorchpi.application.SoftwareComponent, di.unipi.socc.fogtorchpi.infrastructure.ComputationalNode, di.unipi.socc.fogtorchpi.infrastructure.Infrastructure)
	 */
	public double computeNETEnergy(SoftwareComponent s, ComputationalNode n, MobileCloudInfrastructure i) {
		MobileSoftwareComponent cmp = (MobileSoftwareComponent) s;
		MobileCloudInfrastructure mci = (MobileCloudInfrastructure) i;
		return alpha * (cmp.getInData() + cmp.getOutData()) + C; 
	}


	@Override
	public double computeQuantileNETEnergy(Transaction s, ComputationalNode n, MobileCloudInfrastructure i) {
		MobileSoftwareComponent cmp = (MobileSoftwareComponent) s;
		MobileCloudInfrastructure mci = (MobileCloudInfrastructure) i;
		return alpha * (cmp.getInData() + cmp.getOutData()) + C; 
	}

}
