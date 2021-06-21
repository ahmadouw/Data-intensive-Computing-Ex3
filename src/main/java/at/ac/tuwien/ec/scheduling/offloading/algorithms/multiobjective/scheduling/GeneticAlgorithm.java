package at.ac.tuwien.ec.scheduling.offloading.algorithms.multiobjective.scheduling;

import at.ac.tuwien.ec.scheduling.Scheduling;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;

import java.util.ArrayList;

public class GeneticAlgorithm extends OffloadScheduler {
    @Override
    public ArrayList<? extends Scheduling> findScheduling() {
        return null;
    }

    /*
    useful methods to use:

    OffloadScheduler.isValid(OffloadScheduling deployment, MobileSoftwareComponent s, ComputationalNode n)
     checks if MobileSoftwareComponent s can be allocated to node n. It checks if there are enough
     hardware capabilities in the node, if there is enough battery and if connectivity requirements are satisfied

     MobileInfrastructure.getTransmissionTime(MobileSoftwareComponent s, NetworkedNode n, NetworkedNode m):
     computes the actual transmission time to send software component s from node n to node m

     MobileInfrastructure.getDesiredTransmissionTime(MobileSoftwareComponent s, NetworkedNode n, NetworkedNode m):
     computes the transmission time desired by the user, in case there are mobile applications with real-time
     constraints. If yours is not the case, it returns infinite, meaning that the user is happy with
     any transmission time;

    public Set<NetworkConnection> getIncomingLinksTo(ComputationalNode cn):
    returns nodes who have a direct connection with node cn

    MobileSoftwareComponent.getLocalRuntimeOnNode(cn, infrastructur e):
    returns time for local execution (no offloading) on node cn;

    MobileSoftwareComponent.getRuntimeOnNode(n, infrastructure):
    returns time for execution on node n (including time to offload from the mobile device,
    which is identified by the userID of MobileSoftwareComponent.

    OffloadScheduler.deploy (OffloadScheduling deployment, MobileSoftwareComponent s, ComputationalNode n):
    it deploys tasks on computational nodes, performing the necessary calculations (runtime, energy consumption)

    Also, to calculate energy consumption of the nodes, you need to consider two parameters: the CPU
    consumption of the node, i.e., the energy consumption for executing a task locally, and the NET
    consumption of the node, i.e. the energy consumption for transferring the task to another node. Each
    consumption is modelled, respectively, by a CPUEnergyModel object and by a NETEnergyModel
    object. To compute the energy consumption, you need to obtain the required energy model of the
    Mobile device (according to if you are executing locally or remotely) and calculate it using
    CPUEnergyModel.computeCPUEnergy() or NETEnergyModel.computeNETEnergy. Please note: this is
    only the instantaneous consumption, and must be multiplied for the time.

    V: set of nodes (tasks)
    E: set of directed edge (dependencies)
    c: V -> R is a function that associates a weight c(u) to each node u in V
        c(u) represent the execution time of the task Tu which is represented by the node u in V
    tau: function c: E -> R that associates a weight to a directed edge if u and v are two nodes in V
        then tau(u, v) denotes the inter-tasks communication time between Tu and Tv#

    the goal is to minimize the makespan without violation precedence constraints
    where makespan = max{ft(u)}
    st(u) is the start time of task u
    ft(u) is the finish time of task u
     */

    // TODO: 6/21/21 what is my function to optimize is there a priority of the 4 ?

    // todo: generate an initial population

    // todo: implement or call the fitness function for a given cromosom

    // todo: select a pair of individuals based on the fitness function

    // todo: produce next generation from the selected pairs by performing random changes on the slected partners
    //  mutation
    //  crossover

    // todo: stopping criterion
    //  some value ot be satisfied
    //  no improvement for several iterations
    //  better than other algorithm ?

    // todo: decide on which structure to use to represent chromosomes
    // chromosome -> individual in the population
    //  represents a possible solution to a problem
    //  this is a schedule of a group of tasks on a group of processors
    //  can be represented as a sequence of individual schedules -> one for each processor in the group
    //    separated by special values (Tj, Pi)
    //  queue of tasks assigned to that processor

    // todo: inizialise population by randomly placing each task on a processor
    //  what does that exactly mean ?
    //  use different random number generators (Poisson, Nomral, Uniform, Laplace)
    // can lead to configurations that do not correspond to the specification
    // make corrections to make the generated chromosome compliant with the specification or by including penalties in the fitness function

    // todo mutation:
    //  randomly altering certain genes
    //  dynamically adjust the mutation rate depending on the fitness variation
    //   see paper

    // todo: roulett wheel selection method



}
