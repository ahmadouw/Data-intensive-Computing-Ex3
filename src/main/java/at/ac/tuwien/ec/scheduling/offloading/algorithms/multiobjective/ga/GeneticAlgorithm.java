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

    private static final int POPULATION_SIZE = 100;
    private static final int MAX_GENERATIONS = 100;
    private static final double MUTATION_RATE = 0.01;

    private int generation;
    private final IntegerDistribution nodeSampler;

    // metrics
    private final ArrayList<Double> averageFitness = new ArrayList<>();
    private final ArrayList<Double> averageRunTime = new ArrayList<>();
    private final ArrayList<Integer> numberValidSchedules = new ArrayList<>();
    private final ArrayList<Long> fitnessCalculationTime  = new ArrayList<>();

    public GeneticAlgorithm(Tuple2<MobileApplication, MobileCloudInfrastructure> t) {
        super();
        setMobileApplication(t._1());
        setInfrastructure(t._2());
        generation = 0;
        nodeSampler = new UniformIntegerDistribution(0, currentInfrastructure.getAllNodes().size());
    }
    // todo: improve termination criterion
    @Override
    public ArrayList<OffloadScheduling> findScheduling() {

        ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population = generateInitialPopulation();


        while (generation < MAX_GENERATIONS) {

            HashMap<Integer, Double> fitness = calculateFitness(population);
            HashMap<Integer, Double> probabilities = calculateSelectionProbabilities(population, fitness);
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> populationNextGen = new ArrayList<>();
            ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> populationMutated = new ArrayList<>();

            // create new generation with crossover
            while (populationNextGen.size() < POPULATION_SIZE) {
                LinkedHashMap<MobileSoftwareComponent, ComputationalNode> father = rouletteWheelSelection(population, probabilities);
                LinkedHashMap<MobileSoftwareComponent, ComputationalNode> mother = rouletteWheelSelection(population, probabilities);
                if (father != null && mother != null) {
                    List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> children = crossover(father, mother);
                    populationNextGen.addAll(children);
                }
            }

            // apply mutation
            for (LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome : populationNextGen) {
                populationMutated.add(mutate(chromosome));
            }
            population = populationMutated;
            generation++;
        }

        ArrayList<OffloadScheduling> deployments = new ArrayList<>();
        deployments.add(getTopChromosome(population,  calculateFitness(population)));
        return deployments;
    }

    // todo: dynamically adjust the mutation rate, depending on the fitness variation
    private ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> generateInitialPopulation() {
        ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population = new ArrayList<>();

        ArrayList<MobileSoftwareComponent> taskList = getTaskList();
        int initialPopulationSize = POPULATION_SIZE * 10;
        while (population.size() < initialPopulationSize) {

            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome = new LinkedHashMap<>();

            for (MobileSoftwareComponent currTask : taskList) {
                int sample = nodeSampler.sample();
                ComputationalNode node;
                if (currTask.isOffloadable()) {
                    if (sample >= currentInfrastructure.getAllNodes().size()) {
                        node = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
                    } else {
                        node = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[sample];
                    }
                } else {
                    node = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
                }
                chromosome.put(currTask, node);
            }
            population.add(chromosome);
        }
        return population;
    }


    private ArrayList<MobileSoftwareComponent> getTaskList() {
        ArrayList<MobileSoftwareComponent> taskList = new ArrayList<>();
        DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> deps = this.getMobileApplication().getTaskDependencies();
        for (MobileSoftwareComponent dep : deps) {
            taskList.add(dep);
        }
        return taskList;
    }


    private HashMap<Integer, Double> calculateFitness(ArrayList<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> population) {
        HashMap<Integer, Double> fitness = new HashMap<>();
        HashMap<Integer, Double> runTimes = new HashMap<>();

        long startTime = System.currentTimeMillis();

        for (int i = 0; i < population.size(); i++) {
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome = population.get(i);
            OffloadScheduling scheduling = new OffloadScheduling();
            deployAll(scheduling, chromosome);
            if (schedulingValid(scheduling, chromosome)) {
                runTimes.put(i, scheduling.getRunTime());
            } else {
                fitness.put(i, 0.0);
            }
            undeployAll(scheduling, chromosome);
        }

        double minRunTime = Double.MAX_VALUE;
        double maxRunTime = Double.MIN_VALUE;
        for (Map.Entry<Integer, Double> entry : runTimes.entrySet()) {
            double currentRunTime = entry.getValue();
            if (Double.isInfinite(currentRunTime)) {
                continue;
            }
            if (currentRunTime < minRunTime) {
                minRunTime = currentRunTime;
            }
            if (currentRunTime > maxRunTime) {
                maxRunTime = currentRunTime;
            }
        }

        for (Map.Entry<Integer, Double> entry : runTimes.entrySet()) {
            if (!fitness.containsKey(entry.getKey())) {
                double fitnessValue;
                if (Double.isInfinite(entry.getValue())) {
                    fitnessValue = 0.0;
                } else {
                    double f = 1 / entry.getValue();
                    double minFitness = 1 / maxRunTime;
                    double maxFitness = 1 / minRunTime;
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

        long endTime = System.currentTimeMillis();

        // capture metrics
        this.averageFitness.add(calculateAverage(new ArrayList<>(fitness.values())));
        this.averageRunTime.add(calculateAverage(runTimes.values().stream().filter(v -> !Double.isInfinite(v) && v > 0.0).collect(Collectors.toList())));
        this.numberValidSchedules.add((int) fitness.values().stream().filter(v -> v > 0.0).count());
        this.fitnessCalculationTime.add(endTime - startTime);
        return fitness;
    }


    private boolean schedulingValid(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            if (!isValid(scheduling, task, chromosome.get(task))) {
                return false;
            }
        }
        return true;
    }


    private void deployAll(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            deploy(scheduling, task, chromosome.get(task));
        }
    }


    private void undeployAll(OffloadScheduling scheduling, LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {
        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            undeploy(scheduling, task, chromosome.get(task));
        }
    }


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


    private List<LinkedHashMap<MobileSoftwareComponent, ComputationalNode>> crossover(
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> fatherChromosome,
            LinkedHashMap<MobileSoftwareComponent, ComputationalNode> motherChromosome) {

        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> child1 = new LinkedHashMap<>();
        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> child2 = new LinkedHashMap<>();

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

    //  dynamically adjust the mutation rate depending on the fitness variation
    private LinkedHashMap<MobileSoftwareComponent, ComputationalNode> mutate(LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosome) {

        LinkedHashMap<MobileSoftwareComponent, ComputationalNode> chromosomeMutated = new LinkedHashMap<>();

        Set<MobileSoftwareComponent> tasks = chromosome.keySet();
        for (MobileSoftwareComponent task : tasks) {
            ComputationalNode node;
            if (task.isOffloadable()) {
                double mutate = Math.random();
                if (mutate <= MUTATION_RATE) {
                    int sample = nodeSampler.sample();
                    if (sample >= currentInfrastructure.getAllNodes().size()) {
                        node = (ComputationalNode) currentInfrastructure.getNodeById(task.getUserId());
                    } else {
                        node = (ComputationalNode) currentInfrastructure.getAllNodes().toArray()[sample];
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

    private double calculateAverage(List<Double> values) {
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum / values.size();
    }
    
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

}
