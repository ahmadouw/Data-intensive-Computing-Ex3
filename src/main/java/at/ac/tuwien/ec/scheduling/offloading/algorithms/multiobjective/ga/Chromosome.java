package at.ac.tuwien.ec.scheduling.offloading.algorithms.multiobjective.ga;

import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.ComponentLink;
import at.ac.tuwien.ec.model.software.MobileApplication;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduling;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.jgrapht.graph.DirectedAcyclicGraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * The chromosome class is modelled after the RandomOffloadScheduler of the Tutorial
 * Instead of using just a Uniform distribution for creating the scheduling a combination of different
 * random number generators is used as described in the paper Genetic algorithm for DAG scheduling in Grid environments
 * DOI: https://doi.org/10.1109/ICCP.2009.5284747
 */
public class Chromosome extends OffloadScheduler implements Comparable<Chromosome> {

    private int id;
    private OffloadScheduling scheduling;

    // TODO: 6/22/21 its not exactly clear if all nodes are offloaded here or if I should randomly decide to not offload

    public Chromosome(int id, MobileApplication A, MobileCloudInfrastructure I) {
        super();
        this.id = id;
        setMobileApplication(A);
        setInfrastructure(I);
        this.findScheduling();
    }

    public int getId() {
        return id;
    }

    public double getFitness() {
        return scheduling.getRunTime();
    }

    public OffloadScheduling getScheduling() {
        return scheduling;
    }

    @Override
    // override which implements the logic of the scheduling
    public ArrayList<OffloadScheduling> findScheduling() {

        ArrayList<MobileSoftwareComponent> taskList = new ArrayList<>();

        // obtain DAG structure of the workflow modelled in DirectedAcyclicGraph (jgrapht)
        DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> deps = this.getMobileApplication().getTaskDependencies();
        // save nodes in topological order in an auxiliary list
        for (MobileSoftwareComponent dep : deps) {
            taskList.add(dep);
        }

        // in this implementation, the offloading will fail if no offloading target can be found.
        OffloadScheduling scheduling = new OffloadScheduling();
        HashMap<MobileSoftwareComponent, ComputationalNode> targetMap = new HashMap<>();
        ComputationalNode target = null;
        // We implement now the logic of the offloading decisions. First, we iterate over taskList, to
        //select the nodes to be scheduled.
        for(MobileSoftwareComponent currTask : taskList) {
            // Then, if the task is not offloadable, it will try to allocate it
            if(!currTask.isOffloadable()) {
                // on the mobile device which requested execution of the task, which is obtained by using currTask.getUserId().
                // To check if mobile device capabilities are enough for the task, we use the isValid() method.
                if(isValid(scheduling, currTask, (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId())))
                    target = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
            }
            // If currTask can be offloaded
            else {
                // we will offload it to a random computational node.
                // This is done by selecting a random valid target between all computational nodes.
                // todo: use different random number generators (Poisson, Nomral, Uniform, Laplace)
                UniformIntegerDistribution uid = new UniformIntegerDistribution(0, currentInfrastructure.getAllNodes().size() - 1);
                ComputationalNode tmpTarget = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[uid.sample()];
                if(isValid(scheduling, currTask, tmpTarget)) {
                    target = tmpTarget;
                }
            }
            // If a target has been found (if target != null), we deploy the task on the selected target using the deploy method.
            if(target != null) {
                deploy(scheduling, currTask, target);
                targetMap.put(currTask, target);
            }
        }
        // At the end scheduling is added to the schedulings list and returned because findScheduling expects it
        // we always create one random scheduling
        this.scheduling = (OffloadScheduling) scheduling.clone();
        for(MobileSoftwareComponent currTask : taskList) {
            undeploy(scheduling, currTask, targetMap.get(currTask));
        }
        return null;
    }


    @Override
    public int compareTo(Chromosome chromosome) {
        if (scheduling.getRunTime() == chromosome.getScheduling().getRunTime()) {
            return 0;
        } else if (scheduling.getRunTime() < chromosome.getScheduling().getRunTime()) {
            return 1;
        } else {
            return 0;
        }
    }
}
