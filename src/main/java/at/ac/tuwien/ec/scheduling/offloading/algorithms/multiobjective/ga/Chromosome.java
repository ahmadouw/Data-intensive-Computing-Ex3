package at.ac.tuwien.ec.scheduling.offloading.algorithms.multiobjective.ga;

import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.ComponentLink;
import at.ac.tuwien.ec.model.software.MobileApplication;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.scheduling.SchedulingHistogram;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduling;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.jgrapht.graph.DirectedAcyclicGraph;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * The chromosome class is modelled after the RandomOffloadScheduler of the Tutorial
 * Instead of using just a Uniform distribution for creating the scheduling a combination of different
 * random number generators is used as described in the paper Genetic algorithm for DAG scheduling in Grid environments
 * DOI: https://doi.org/10.1109/ICCP.2009.5284747
 */
public class Chromosome extends OffloadScheduler {

    // TODO: 6/22/21 its not exactly clear if all nodes are offloaded here or if I should randomly decide to not offload

    public Chromosome(MobileApplication A, MobileCloudInfrastructure I) {
        super();
        setMobileApplication(A);
        setInfrastructure(I);
    }


    // same as first constructor but input within a tuple
    public Chromosome(Tuple2<MobileApplication,MobileCloudInfrastructure> t) {
        super();
        setMobileApplication(t._1());
        setInfrastructure(t._2());
    }

    @Override
    // override which implements the logic of the scheduling
    public ArrayList<OffloadScheduling> findScheduling() {
        // declare empty list of scheduling which is the output of the method
        ArrayList<OffloadScheduling> schedulings = new ArrayList<OffloadScheduling>();

        ArrayList<MobileSoftwareComponent> taskList = new ArrayList<MobileSoftwareComponent>();

        // obtain DAG structure of the workflow modelled in DirectedAcyclicGraph (jgrapht)
        DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> deps = this.getMobileApplication().getTaskDependencies();
        // save nodes in topological order in an auxiliary list
        Iterator<MobileSoftwareComponent> it = deps.iterator();
        while(it.hasNext()) {
            taskList.add(it.next());
        }

        // in this implementation, the offloading will fail if no offloading target can be found.
        OffloadScheduling scheduling = new OffloadScheduling();
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
                UniformIntegerDistribution uid = new UniformIntegerDistribution(0, currentInfrastructure.getAllNodes().size() - 1);
                ComputationalNode tmpTarget = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[uid.sample()];
                if(isValid(scheduling, currTask, tmpTarget)) {
                    target = tmpTarget;
                }
            }
            // If a target has been found (if target != null), we deploy the task on the selected target using the deploy method.
            if(target != null) {
                deploy(scheduling, currTask, target);
            }
            // At the end of the loop, scheduling is added to the schedulings list and returned as output.
            schedulings.add(scheduling);
        }
            return schedulings;
    }



}
