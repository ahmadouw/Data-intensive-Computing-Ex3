package at.ac.tuwien.ec.model.infrastructure.energy;

import java.io.Serializable;

import at.ac.tuwien.ec.blockchain.Transaction;
import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.model.software.SoftwareComponent;

public class MobileWiFiNETEnergyModel implements NETEnergyModel,Serializable {

	private final double alpha = 0.007e-3;
	private final double C = 5.9e-3;
	
	/* We assume that inData and outData are in bytes
	 * (non-Javadoc)
	 * @see at.ac.tuwien.ec.infrastructuremodel.energy.NETEnergyModel#computeNETEnergy(di.unipi.socc.fogtorchpi.application.SoftwareComponent, di.unipi.socc.fogtorchpi.infrastructure.ComputationalNode, di.unipi.socc.fogtorchpi.infrastructure.Infrastructure)
	 */
	public double computeNETEnergy(SoftwareComponent s, ComputationalNode n, MobileCloudInfrastructure i) {
		MobileSoftwareComponent cmp = (MobileSoftwareComponent) s;
		MobileCloudInfrastructure mci = (MobileCloudInfrastructure) i;
		ComputationalNode mDevice = (ComputationalNode) i.getNodeById(s.getUserId());
		double utilization = (cmp.getOutData())/mci.getLinkBandwidth(mDevice,n);
		return (alpha * utilization + C) * mci.getTransmissionTime((MobileSoftwareComponent) s,mDevice,n) ; 
	}

	@Override
	public double computeQuantileNETEnergy(Transaction s, ComputationalNode n, MobileCloudInfrastructure i) {
		MobileSoftwareComponent cmp = (MobileSoftwareComponent) s;
		MobileCloudInfrastructure mci = (MobileCloudInfrastructure) i;
		ComputationalNode mDevice = (ComputationalNode) i.getNodeById(s.getUserId());
		double utilization = (cmp.getOutData())/mci.getLinkBandwidth(mDevice,n);
		return (alpha * utilization + C) * mci.getTransmissionTime((MobileSoftwareComponent) s,mDevice,n) ; 
	}

}
