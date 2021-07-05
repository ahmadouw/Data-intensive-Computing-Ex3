package at.ac.tuwien.ec.scheduling.offloading.algorithms.heuristics;


import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.jgrapht.Graph;
import org.jgrapht.graph.DirectedAcyclicGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;

import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.MobileDevice;
import at.ac.tuwien.ec.model.software.ComponentLink;
import at.ac.tuwien.ec.model.software.MobileApplication;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduling;
import at.ac.tuwien.ec.scheduling.offloading.algorithms.heftbased.utils.NodeRankComparator;
import at.ac.tuwien.ec.scheduling.utils.RuntimeComparator;
import at.ac.tuwien.ec.sleipnir.OffloadingSetup;
import scala.Tuple2;
import sun.security.provider.certpath.Vertex;

/**
 * OffloadScheduler class that implements the
 * Heterogeneous Earliest-Finish-Time (HEFT) algorithm
 * , a static scheduling heuristic, for efficient application scheduling
 *
 * H. Topcuoglu, S. Hariri and Min-You Wu,
 * "Performance-effective and low-complexity task scheduling for heterogeneous computing,"
 * in IEEE Transactions on Parallel and Distributed Systems, vol. 13, no. 3, pp. 260-274, March 2002, doi: 10.1109/71.993206.
 */

public class HLFET extends OffloadScheduler {
    /**
     *
     * @param A MobileApplication property from  SimIteration
     * @param I MobileCloudInfrastructure property from  SimIteration
     * Constructors set the parameters and calls setRank() to nodes' ranks
     */

    public HLFET(MobileApplication A, MobileCloudInfrastructure I) {
        super();
        setMobileApplication(A);
        setInfrastructure(I);
        setRank(this.currentApp,this.currentInfrastructure);
    }

    public HLFET(Tuple2<MobileApplication,MobileCloudInfrastructure> t) {
        super();
        setMobileApplication(t._1());
        setInfrastructure(t._2());
        setRank(this.currentApp,this.currentInfrastructure);
    }

    /**
     * Processor selection phase:
     * select the tasks in order of their priorities and schedule them on its "best" processor,
     * which minimizes task's finish time
     * @return
     */
    @Override
    public ArrayList<? extends OffloadScheduling> findScheduling() {
        double start = System.nanoTime();
        /*scheduledNodes contains the nodes that have been scheduled for execution.
         * Once nodes are scheduled, they are taken from the PriorityQueue according to their runtime
         */
        PriorityQueue<MobileSoftwareComponent> scheduledNodes
                = new PriorityQueue<MobileSoftwareComponent>(new RuntimeComparator());
        /*
         * tasks contains tasks that have to be scheduled for execution.
         * Tasks are selected according to their upRank (at least in HEFT)
         */
        PriorityQueue<MobileSoftwareComponent> tasks = new PriorityQueue<MobileSoftwareComponent>(new NodeRankComparator());
        //To start, we add all nodes in the workflow
        tasks.addAll(currentApp.getTaskDependencies().vertexSet());
        ArrayList<OffloadScheduling> deployments = new ArrayList<OffloadScheduling>();

        MobileSoftwareComponent currTask;
        //We initialize a new OffloadScheduling object, modelling the scheduling computer with this algorithm
        OffloadScheduling scheduling = new OffloadScheduling();
        //We check until there are nodes available for scheduling
        while((currTask = tasks.peek()) != null)
        {
            double tMin = Double.MAX_VALUE; //Earliest starting time for the task
            ComputationalNode target = null;

            if(!currTask.isOffloadable())
            {
                // If task is not offloadable, deploy it in the mobile device (if enough resources are available)
                if(isValid(scheduling,currTask,(ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId())))
                    target = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());

            }
            else
            {
                //Check the earliest starting time for all available Cloud/Edge nodes
                for(ComputationalNode cn : currentInfrastructure.getAllNodes())
                    if(cn.getESTforTask(currTask) < tMin &&
                            isValid(scheduling,currTask,cn))
                    {
                        tMin = cn.getESTforTask(currTask) ; // Earliest Starting Time EST;
                        target = cn;

                    }
                /*
                 * We need this check, because there are cases where, even if the task is offloadable,
                 * local execution is the best option
                 */
                ComputationalNode localDevice = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
                if(localDevice.getESTforTask(currTask) < tMin &&
                        isValid(scheduling,currTask,localDevice))
                {
                    target = localDevice;

                }

            }
            //if scheduling found a target node for the task, it allocates it to the target node
            if(target != null)
            {
                deploy(scheduling,currTask,target);
                scheduledNodes.add(currTask);
                tasks.remove(currTask);
            }
            else if(!scheduledNodes.isEmpty());
            {
                MobileSoftwareComponent terminated = scheduledNodes.remove();
                ((ComputationalNode) scheduling.get(terminated)).undeploy(terminated);
            }
            /*
             * if simulation considers mobility, perform post-scheduling operations
             * (default is to update coordinates of mobile devices)
             */
            if(OffloadingSetup.mobility)
                postTaskScheduling(scheduling);
        }
        double end = System.nanoTime();
        scheduling.setExecutionTime(end-start);
        deployments.add(scheduling);
        return deployments;
    }

    protected void setRank(MobileApplication A, MobileCloudInfrastructure I)
    {
        //for the calculation of b-levels, the reverse topological order is needed
        Set<MobileSoftwareComponent> topological = A.getTaskDependencies().vertexSet();
        List<MobileSoftwareComponent> reversed = new ArrayList<>(topological);
        Collections.reverse(reversed);
        for(MobileSoftwareComponent msc : reversed)
            bLevel(msc,A.getTaskDependencies(),I);

    }

    /**
     * bLevel is the task prioritizing phase of HLFET
     * @param msc
     * @param dag Mobile Application's DAG
     * @param infrastructure
     * @return the static b-level of msc
     */
    private double bLevel(MobileSoftwareComponent msc, DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> dag,
                          MobileCloudInfrastructure infrastructure) {
        double w_cmp = 0.0; // average execution time of task on each processor / node of this component, which is used as the cost of the current task
            int numberOfNodes = infrastructure.getAllNodes().size() + 1;
            for(ComputationalNode cn : infrastructure.getAllNodes())
                w_cmp += msc.getLocalRuntimeOnNode(cn, infrastructure);

            w_cmp = w_cmp / numberOfNodes;

            double tmpBLevel;
            double max = 0; // max of child b-level + edge weight from current task to child
            // as we start at the last node, make sure that the b-level is just it's average execution time (as it does not have any children)
            if(dag.outDegreeOf(msc) == 0){
                msc.setRank(w_cmp);
            }else {
                for (ComponentLink neigh : dag.outgoingEdgesOf(msc))
                {
                    // b-level = w(ni) +  max(c(ni,ny)+ b-level(ny)
                    // where c(ni,ny) is the average commmunication cost of edge (i, y)
                    tmpBLevel = neigh.getTarget().getRank(); // child's b-level
                    double tmpCost = 0;  //
                    // we consider only offloadable children. If a child is not offloadable, communication cost is 0
                    if (neigh.getTarget().isOffloadable()) {
                        // get the average communication cost
                        for (ComputationalNode cn : infrastructure.getAllNodes())
                            tmpCost += infrastructure.getTransmissionTime(neigh.getTarget(), infrastructure.getNodeById(msc.getUserId()), cn);
                        tmpCost = tmpCost / (infrastructure.getAllNodes().size());
                    }
                    // calculate the sum of child's b-level and the average cost to the childaxy
                    double tmpRank = tmpBLevel + tmpCost;
                    max = (tmpRank > max) ? tmpRank : max;
                }
                msc.setRank(w_cmp + max);
            }

        return msc.getRank();
    }

}
