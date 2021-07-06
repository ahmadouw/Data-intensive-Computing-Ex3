package at.ac.tuwien.ec.scheduling.offloading.algorithms.multiobjective.ga;

import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.ComponentLink;
import at.ac.tuwien.ec.model.software.MobileApplication;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduling;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.jgrapht.graph.DirectedAcyclicGraph;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

public class GeneticAlgorithm extends OffloadScheduler {

    // Hyperparameters of the GA
    private static final int POPULATION_SIZE = 100;
    private static final int MAX_GENERATIONS = 100;
    private static final double MUTATION_RATE = 0.01;
    private static final double ELITISM_RATE = 0.05;
    private static final double RATE_LOCAL_COMPUTATION = 0.2;
    private static final double WEIGHT_RUN_TIME = 0.5;
    private static final double WEIGHT_BATTERY_LIFE = 0.5;

    private int generation;
    private final IntegerDistribution nodeSampler;

    // metrics captured over generations
    private final ArrayList<Double> averageFitness = new ArrayList<>();
    private final ArrayList<Double> averageRunTime = new ArrayList<>();
    private final ArrayList<Double> averageBatteryLife = new ArrayList<>();
    private final ArrayList<Integer> numberValidSchedules = new ArrayList<>();
    private final ArrayList<Long> fitnessCalculationTime  = new ArrayList<>();

    public GeneticAlgorithm(Tuple2<MobileApplication, MobileCloudInfrastructure> t) {
        super();
        setMobileApplication(t._1());
        setInfrastructure(t._2());
        generation = 0;
        nodeSampler = new UniformIntegerDistribution(0, currentInfrastructure.getAllNodes().size() - 1);
    }


    @Override
    public ArrayList<OffloadScheduling> findScheduling() {

        // generate initial population
        ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population = generateInitialPopulation();

        // optimizer over multiple generations
        while (generation < MAX_GENERATIONS) {
            // next population
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> populationNextGen = new ArrayList<>();

            // calculate fitness values for current population
            HashMap<Integer, Double> fitness = calculateFitness(population);

            // elitism
            populationNextGen.addAll(getElites(population, fitness));

            // get probabilities for chromosomes for being selected for crossover
            HashMap<Integer, Double> probabilities = calculateSelectionProbabilities(population, fitness);
            // create new chromosomes with crossover until POPULATION_SIZE reached
            while (populationNextGen.size() < POPULATION_SIZE) {
                // select parents using roulette wheel selection
                LinkedHashMap<MobileSoftwareComponent, ComputationalNode> father = rouletteWheelSelection(population, probabilities);
                LinkedHashMap<MobileSoftwareComponent, ComputationalNode> mother = rouletteWheelSelection(population, probabilities);
                if (father != null && mother != null) {
                    // create children
                    List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> children = crossover(father, mother);
                    populationNextGen.addAll(children);
                }
            }

            // apply mutation to next generation
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> populationMutated = new ArrayList<>();
            for (LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome : populationNextGen) {
                populationMutated.add(mutate(chromosome));
            }
            population = populationMutated;
            generation++;
        }

        // return the chromosome with the top fitness of the last population
        ArrayList<OffloadScheduling> deployments = new ArrayList<>();
        deployments.add(getTopChromosome(population,  calculateFitness(population)));
        return deployments;
    }

   /*
   Generates the initial population randomly using Uniform distribution.
   This is modelled after the RandomOffloading example provided.
    */
    private ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> generateInitialPopulation() {
        ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population = new ArrayList<>();

        // get task list to ensure DAG constraints ensured
        ArrayList<MobileSoftwareComponent> taskList = getTaskList();

        // the initial population size is larger than the subsequent populations
        // this is because a not valid scheduling can be created randomly and only valid schedules are considered during selection
        int initialPopulationSize = POPULATION_SIZE * 10;
        while (population.size() < initialPopulationSize) {
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome = new LinkedHashMap<>();

            // iterate all tasks in correct order
            for (MobileSoftwareComponent currTask : taskList) {
                ComputationalNode node;
                // if task can be offloaded we randomly decide if task is computed locally or offloaded
                // then the node is randomly selected from all cloud and edge nodes
                if (currTask.isOffloadable()) {
                    double offload = Math.random();
                    if (offload > RATE_LOCAL_COMPUTATION) {
                        node = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[nodeSampler.sample()];
                    } else {
                        node = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
                    }
                }
                // if task can not be offloaded we schedule it to the correct mobile device
                else {
                    node = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
                }
                chromosome.put(currTask, node);
            }
            population.add(chromosome);
        }
        return population;
    }


    /*
    Create the task list from the DAG to ensure tasks are scheduled in correct order
     */
    private ArrayList<MobileSoftwareComponent> getTaskList() {
        ArrayList<MobileSoftwareComponent> taskList = new ArrayList<>();
        DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> deps = this.getMobileApplication().getTaskDependencies();
        for (MobileSoftwareComponent dep : deps) {
            taskList.add(dep);
        }
        return taskList;
    }


    /*
    Fitness function:
    We optimize for run time and battery life of a scheduling. (mentioned in Forum)
    It is necessary to schedule all tasks to calculate the fitness function.
    The fitnessFunction can be weighted to prioritize one of these values. (weights must sum to 1.0)
    For both values we normalize to range [0, 1].
    For run time the inverse is computed to measure minimizing runtime.
    For battery life we try to maximize (not inverted).
     */
    private HashMap<Integer, Double> calculateFitness(ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population) {
        HashMap<Integer, Double> fitness = new HashMap<>();
        HashMap<Integer, Double> notValid = new HashMap<>();
        HashMap<Integer, Double> runTimes = new HashMap<>();
        HashMap<Integer, Double> batteryLife = new HashMap<>();

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < population.size(); i++) {
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome = population.get(i);
            OffloadScheduling scheduling = new OffloadScheduling();
            deployAll(scheduling, chromosome);
            if (schedulingValid(scheduling, chromosome)) {
                runTimes.put(i, scheduling.getRunTime());
                batteryLife.put(i, scheduling.getBatteryLifetime());
            } else {
                notValid.put(i, 0.0);
            }
            undeployAll(scheduling, chromosome);
        }

        HashMap<Integer, Double> fitnessRunTime = fitness(notValid, runTimes, true);
        HashMap<Integer, Double> fitnessBatteryLife = fitness(notValid, batteryLife, false);
        for (int i = 0; i < population.size(); i++) {
            if (notValid.containsKey(i)) {
                fitness.put(i, 0.0);
            } else {
                fitness.put(i, (fitnessRunTime.get(i) * WEIGHT_RUN_TIME) + (fitnessBatteryLife.get(i) * WEIGHT_BATTERY_LIFE));
            }
        }

        long endTime = System.currentTimeMillis();

        // capture metrics
        this.averageFitness.add(calculateAverage(new ArrayList<>(fitness.values())));
        this.averageRunTime.add(calculateAverage(runTimes.values().stream().filter(v -> !Double.isInfinite(v) && v > 0.0).collect(Collectors.toList())));
        this.averageBatteryLife.add(calculateAverage(batteryLife.values().stream().filter(v -> !Double.isInfinite(v) && v > 0.0).collect(Collectors.toList())));
        this.numberValidSchedules.add((int) fitness.values().stream().filter(v -> v > 0.0).count());
        this.fitnessCalculationTime.add(endTime - startTime);
        return fitness;
    }


    /*
    Helper function to compute a normalized value in range [0, 1] for the values to optimize.
     */
    private HashMap<Integer, Double> fitness(HashMap<Integer, Double> notValid, HashMap<Integer, Double> map, boolean inverse) {
        HashMap<Integer, Double> fitness = new HashMap<>();

        ArrayList<Double> values = new ArrayList<>(map.values());
        double min = min(values);
        double max = max(values);

        for (Map.Entry<Integer, Double> entry : map.entrySet()) {
            if (notValid.containsKey(entry.getKey())) {
                fitness.put(entry.getKey(), notValid.get(entry.getKey()));
            } else {
                double fitnessValue;
                if (Double.isInfinite(entry.getValue())) {
                    fitnessValue = 0.0;
                } else {
                    double f, minFitness, maxFitness;
                    if (inverse) {
                        f = 1 / entry.getValue();
                        minFitness = 1 / max;
                        maxFitness = 1 / min;
                    } else {
                        f = entry.getValue();
                        minFitness = min;
                        maxFitness = max;
                    }
                    double normalized = (f - minFitness) / (maxFitness - minFitness);
                    if (Double.isNaN(normalized)) {
                        fitnessValue = 0.0;
                    } else {
                        fitnessValue = normalized;
                    }
                }
                fitness.put(entry.getKey(), fitnessValue);
            }
        }
        return fitness;
    }


    /*
    Elitism:
    Function to select the top chromosomes by their fitness.
    Elites survive as they are and continue to the next generation without crossover.
     */
    private List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> getElites(
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population,
            HashMap<Integer, Double> fitness) {

        List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> elites = new ArrayList<>();
        int numElites = (int) (POPULATION_SIZE * ELITISM_RATE);

        List<Double> fitnessValues = fitness.values().stream().sorted().collect(Collectors.toList());
        double fitnessCutoff = fitnessValues.get(fitness.size() - numElites);

        for (int i = 0; i < population.size(); i++) {
            if (fitness.get(i) >= fitnessCutoff) {
                elites.add(population.get(i));
            }
        }

        return elites;
    }


    /*
    Selection function:
    Implementation of the roulette wheel selection to randomly pick chromosomes for the next generation.
    In the roulette wheel selection the probability of being selected is proportional to the fitness of the chromosome.
    Not valid schedules (fitness == 0.0) have no chance of being selected.
     */
    private LinkedHashMap<MobileSoftwareComponent, ComputationalNode> rouletteWheelSelection(
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population,
            HashMap<Integer, Double> probabilities) {

        double partialSum = 0.0;
        double roulette = Math.random();

        for (int i = 0; i < population.size(); i++) {
            partialSum += probabilities.get(i);
            if (partialSum >= roulette) {
                return population.get(i);
            }
        }
        return null;
    }


    /*
    Crossover function:
    Two parent chromosomes are used to pass over there genes to the next generation and create two child chromosomes.
    One-point crossover is used for this implementation.
     */
    private List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> crossover(
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> fatherChromosome,
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> motherChromosome) {

        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> child1 = new LinkedHashMap<>();
        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> child2 = new LinkedHashMap<>();

        // random point in chromosome where crossover happens
        IntegerDistribution crossoverPointSampler = new UniformIntegerDistribution(0, fatherChromosome.size() - 1);
        int sample = crossoverPointSampler.sample();

        int taskNumber = 0;
        Set<MobileSoftwareComponent> tasks = fatherChromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            if (taskNumber < sample) {
                child1.put(task, fatherChromosome.get(task));
                child2.put(task, motherChromosome.get(task));
            } else {
                child1.put(task, motherChromosome.get(task));
                child2.put(task, fatherChromosome.get(task));
            }
            taskNumber++;
        }
        return Arrays.asList(child1, child2);
    }


    /*
    Mutation function.
    Randomly mutates chromosomes of the new generation before they are passed on to the next generation.
    Each task in the scheduling has a small chance to be scheduled to another node.
    The approach of randomly scheduling is similar to the one used for the initial creation of the population.
    Mutation can create a scheduling that is not valid!
     */
    private LinkedHashMap<MobileSoftwareComponent, ComputationalNode> mutate(LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {

        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosomeMutated = new LinkedHashMap<>();

        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            ComputationalNode node;
            if (task.isOffloadable()) {
                double mutate = Math.random();
                if (mutate <= MUTATION_RATE) {
                    double offload = Math.random();
                    if (offload > RATE_LOCAL_COMPUTATION) {
                        node = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[nodeSampler.sample()];
                    } else {
                        node = (ComputationalNode) currentInfrastructure.getNodeById(task.getUserId());
                    }
                } else {
                    node = chromosome.get(task);
                }
            } else {
                node = chromosome.get(task);
            }
            chromosomeMutated.put(task, node);
        }
        return chromosomeMutated;
    }


    /*
    Check if a created scheduling is valid.
     */
    private boolean schedulingValid(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            if (!isValid(scheduling, task, chromosome.get(task))) {
                return false;
            }
        }
        return true;
    }


    /*
    Deploy all tasks in the chromosome to the OffloadScheduling.
     */
    private void deployAll(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            deploy(scheduling, task, chromosome.get(task));
        }
    }


    /*
    Undeploy all tasks in the chromosome from the OffloadScheduling to release resources.
    This is necessary because the Infrastructure is shared over the whole population.
    Thus, we have to undeploy after the fitness calculation of each Chromosome.
     */
    private void undeployAll(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            undeploy(scheduling, task, chromosome.get(task));
        }
    }


    /*
    Calculates the probabilities of being selected by the roulette wheel selection.
    The fitness defines the probability of being selected.
     */
    private HashMap<Integer, Double> calculateSelectionProbabilities(
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population,
            HashMap<Integer, Double> fitness) {

        HashMap<Integer, Double> probabilities = new HashMap<>();

        double totalSum = 0.0;
        for (int i = 0; i < population.size(); i++) {
            totalSum += fitness.get(i);
        }

        for (int i = 0; i < population.size(); i++) {
            if (fitness.get(i) > 0.0) {
                probabilities.put(i, fitness.get(i) / totalSum);
            } else {
                probabilities.put(i, 0.0);
            }
        }

        return probabilities;
    }


    /*
    Create an OffloadScheduling from the top chromosome (highest fitness) and deploy all tasks to infrastructure before returning it.
     */
    private OffloadScheduling getTopChromosome(ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population, HashMap<Integer, Double> fitness) {
        OffloadScheduling scheduling = new OffloadScheduling();
        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> topChromosome = new LinkedHashMap<>();

        double maxFitness = Double.MIN_VALUE;
        for (int i = 0; i < population.size(); i++) {
            if (fitness.get(i) > maxFitness) {
                maxFitness = fitness.get(i);
                topChromosome = population.get(i);
            }
        }

        deployAll(scheduling, topChromosome);
        return scheduling;
    }


    /*
    Util function to get minimum value
     */
    private double min(ArrayList<Double> values) {
        double min = Double.MAX_VALUE;

        for (double value : values) {
            if (Double.isInfinite(value)) {
                continue;
            }
            if (value < min) {
                min = value;
            }
        }
        return min;
    }


    /*
    Util function to get maximum value
     */
    private double max(ArrayList<Double> values) {
        double max = Double.MIN_VALUE;

        for (double value : values) {
            if (Double.isInfinite(value)) {
                continue;
            }
            if (value > max) {
                max = value;
            }
        }
        return max;
    }


    /*
    Util function to calculate the average
     */
    private double calculateAverage(List<Double> values) {
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum / values.size();
    }

}
