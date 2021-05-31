package at.ac.tuwien.ec.sleipnir.fgcs;

import org.apache.commons.lang.math.RandomUtils;

import at.ac.tuwien.ec.model.Hardware;
import at.ac.tuwien.ec.model.HardwareCapabilities;
import at.ac.tuwien.ec.model.infrastructure.energy.AMDCPUEnergyModel;
import at.ac.tuwien.ec.model.infrastructure.energy.CPUEnergyModel;
import at.ac.tuwien.ec.model.infrastructure.energy.Mobile3GNETEnergyModel;
import at.ac.tuwien.ec.model.infrastructure.energy.NETEnergyModel;
import at.ac.tuwien.ec.model.infrastructure.energy.SamsungS2DualEnergyModel;
import at.ac.tuwien.ec.model.pricing.EdgePricingModel;

public class FGCSSetup {
	
	public static void readFromFile(String url){
		
	} 
	
	public static double workflowMips = 20000;
	public static double workflowIndata = 250000;
	public static double workflowOutData = 250000;
	
	public static RandomUtils rand = new RandomUtils();
	public static int MAP_M = 3;
	public static int MAP_N = 2;
	public static int cloudMaxHops = 1;
	public static int cloudCoreNum = 64;
	public static double cloudRam = 256;
	public static double cloudStorage = 5000;
	public static double cloudMipsPerCore = 100000.0;
	public static int edgeCoreNum = 16;
	public static double edgeRam = 128;
	public static double edgeStorage = 1000;
	public static double edgeMipsPerCore = 80000.0;
	public static EdgePricingModel edgePricingModel = new EdgePricingModel();
	public static CPUEnergyModel edgeCPUEnergyModel = new AMDCPUEnergyModel();
	public static int mobileNum = 2;
	public static double mobileEnergyBudget = 10000;
	public static HardwareCapabilities defaultMobileDeviceHardwareCapabilities = 
			new HardwareCapabilities(new Hardware(2,16,(int)16e10),0.0);
	public static CPUEnergyModel defaultMobileDeviceCPUModel = new SamsungS2DualEnergyModel();
	public static NETEnergyModel defaultMobileDeviceNETModel = new Mobile3GNETEnergyModel();
	public static double wifiAvailableProbability = 1.0;
	public static HardwareCapabilities defaultEdgeNodeCapabilities = new HardwareCapabilities(
			new Hardware(FGCSSetup.edgeCoreNum,
			FGCSSetup.edgeRam,
			FGCSSetup.edgeStorage),
			FGCSSetup.edgeMipsPerCore);
	public static HardwareCapabilities defaultCloudCapabilities = new HardwareCapabilities(
			new Hardware(FGCSSetup.cloudCoreNum,
			FGCSSetup.cloudRam,
			FGCSSetup.cloudStorage),
			FGCSSetup.cloudMipsPerCore);
	public static double antivirusFileSize;
	public static double task_multiplier = 1;
	public static int lambdaLatency = 0;
	public static int chessMovesNum = 1;
	public static double chess_mi = 10e3;
	public static double facebookImageSize = 20e3;
	public static double facerecImageSize = 10e3;
	public static double navigatorMapSize = 25e6;
	public static String[] algorithms;
	public static boolean batch;
	public static double batteryCapacity = mobileEnergyBudget * mobileNum;
	public static int iterations = 1000;
	public static int cloudNum = 6;
	public static double EchoGamma;
	public static double EchoAlpha;
	public static double EchoBeta;
	public static double Eta = 1.0;
	public static int appNumber = 25;
	

}
