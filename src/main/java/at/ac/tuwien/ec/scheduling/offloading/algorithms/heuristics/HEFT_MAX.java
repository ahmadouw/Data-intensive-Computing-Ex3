package at.ac.tuwien.ec.scheduling.algorithms.heftbased;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.jgrapht.Graph;
import org.jgrapht.graph.DirectedAcyclicGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;

import at.ac.tuwien.ec.model.infrastructure.MobileCloudInfrastructure;
import at.ac.tuwien.ec.model.infrastructure.computationalnodes.ComputationalNode;
import at.ac.tuwien.ec.model.software.ComponentLink;
import at.ac.tuwien.ec.model.software.MobileApplication;
import at.ac.tuwien.ec.model.software.MobileSoftwareComponent;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduling;
import at.ac.tuwien.ec.scheduling.offloading.OffloadScheduler;
import at.ac.tuwien.ec.scheduling.offloading.algorithms.heftbased.utils.NodeRankComparator;
import at.ac.tuwien.ec.scheduling.utils.RuntimeComparator;
import scala.Tuple2;


public class HEFT_MAX extends OffloadScheduler { //extends the OffloadScheduler
	
	public HEFT_MAX(MobileApplication A, MobileCloudInfrastructure I) {
		super(); //call the superclass methods and access the superclass constructor
		setMobileApplication(A);
		setInfrastructure(I);
		setRank(this.currentApp,this.currentInfrastructure); //function will be defined later 
	}
	
	public HEFT_MAX(Tuple2<MobileApplication,MobileCloudInfrastructure> t) { 
		super();
		setMobileApplication(t._1());
		setInfrastructure(t._2());
		setRank(this.currentApp,this.currentInfrastructure); //function will be defined later 
	}
	
	
	@Override //needs to override findscheduling() method (which implements the scheduling logic)
	public ArrayList<OffloadScheduling> findScheduling() {
		// output is an ArrayList of OFFLOADSCHEDULING (needs to stay this way!)
		double start = System.nanoTime();
		PriorityQueue<MobileSoftwareComponent> scheduledNodes //new priorityqueue called "scheduleNodes" of type MobileSoftwareComponent
		= new PriorityQueue<MobileSoftwareComponent>(new RuntimeComparator());
		//ArrayList<MobileSoftwareComponent> tasks = new ArrayList<MobileSoftwareComponent>();
		PriorityQueue<MobileSoftwareComponent> tasks = new PriorityQueue<MobileSoftwareComponent>(new NodeRankComparator());
		//new priorityqueue called "tasks" of type MobileSoftwareComponent
		//important: priority queue is filled with NodeRankComparator, i.e. with the sorted ranks 
		tasks.addAll(currentApp.getTaskDependencies().vertexSet()); //add all vertex (=nodes)= tasks to the tasks
		//Collections.sort(tasks, new NodeRankComparator());
		ArrayList<OffloadScheduling> deployments = new ArrayList<OffloadScheduling>(); 

		tasks.addAll(currentApp.getTaskDependencies().vertexSet()); //add all of the nodes a second time
		
		double currentRuntime;
		MobileSoftwareComponent currTask; 
		OffloadScheduling scheduling = new OffloadScheduling(); 
		while((currTask = tasks.poll())!=null) //until all tasks are queued (task queue is empty) 
		{ 
			//currTask = tasks.remove(0);
			
			if(!scheduledNodes.isEmpty())//if the scheduled nodes are not empty (i.e. some values are in there)
			{
				MobileSoftwareComponent firstTaskToTerminate = scheduledNodes.remove(); //will remove first element of priority queue (scheduled nodes)
				currentRuntime = firstTaskToTerminate.getRunTime(); //get the run time from this software component
				//currentApp.removeEdgesFrom(firstTaskToTerminate);
				//currentApp.removeTask(firstTaskToTerminate);
				((ComputationalNode) scheduling.get(firstTaskToTerminate)).undeploy(firstTaskToTerminate); //undeploy this certain task from the scheduling list
				//scheduledNodes.remove(firstTaskToTerminate);
			}
			double tMin = Double.MAX_VALUE; 
			ComputationalNode target = null; 
			if(!currTask.isOffloadable())  //if this task is not offloadable
			{
				if(isValid(scheduling,currTask,(ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId())))
				{
					target = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
					scheduledNodes.add(currTask); //add the  task to the scheduled nodes (priority queue)
				}
				else
				{
					if(scheduledNodes.isEmpty())
						target = null;
				}
			}
			else //task can be offloaded (potentially!) 
			{
				double maxP = Double.MIN_VALUE; 
				for(MobileSoftwareComponent cmp : currentApp.getPredecessors(currTask)) //loop through all predecessors of the current task
					if(cmp.getRunTime()>maxP)  
						maxP = cmp.getRunTime(); //raise threshold "maxP"
				
				for(ComputationalNode cn : currentInfrastructure.getAllNodes()) //loop through ALL tasks 
					if(maxP + currTask.getRuntimeOnNode(cn, currentInfrastructure) < tMin &&
							isValid(scheduling,currTask,cn))  
					{
						tMin = maxP + currTask.getRuntimeOnNode(cn, currentInfrastructure); //define new threshold
						target = cn; 
					}
				if(maxP + currTask.getRuntimeOnNode((ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId()), currentInfrastructure) < tMin
						&& isValid(scheduling,currTask,(ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId())))
					target = (ComputationalNode) currentInfrastructure.getNodeById(currTask.getUserId());
				// task with lowest threshold is defined as "target" at the end 				
			}
			if(target != null) //there still is another potential target task
			{
				deploy(scheduling,currTask,target); //deploy current task
				scheduledNodes.add(currTask); // current task is added to the scheduled nodes again 
			}
			else //no potential target task 
			{
				if(scheduledNodes.isEmpty()) // no scheduled tasks
					target = null; 
			}
								
		}
		double end = System.nanoTime(); // end time 
		scheduling.setExecutionTime(end-start); //calculate execution time 
		deployments.add(scheduling); // add scheduling to deployment 
		return deployments; // return ArrayList<OffloadScheduling> named deployments
	}

	private void setRank(MobileApplication A, MobileCloudInfrastructure I) //define function for setRank 
	{
		for(MobileSoftwareComponent msc : A.getTaskDependencies().vertexSet()) //get all tasks 
			msc.setVisited(false); // set every task to "not yet visited"
				
		for(MobileSoftwareComponent msc : A.getTaskDependencies().vertexSet()) //same
			upRank(msc,A.getTaskDependencies(),I); //calculate upward rank
				
	}

	//calculate upward-rank =>  important to assign a priority to every task (!)
	private double upRank(MobileSoftwareComponent msc, DirectedAcyclicGraph<MobileSoftwareComponent, ComponentLink> dag,
			MobileCloudInfrastructure I) { 
			
		if(!msc.isVisited()) //if component has not yet been visited
		{
			msc.setVisited(true); // define it as visited  
			int numberOfNodes = I.getAllNodes().size() + 1; // number of nodes + 1 Anzahl der Nodes in der Infrastruktur + 1 (für running time auf aktueller Node)

			List<Double> f1 = new ArrayList<Double>(); //new arraylist to store the ranks 
			for(ComputationalNode cn : I.getAllNodes()) { //loop through all nodes 
				f1.add(msc.getLocalRuntimeOnNode(cn, I)); //get runtime of task on every computational node in infrastructure and add it to f1
			}

			f1.add(msc.getLocalRuntimeOnNode((ComputationalNode) I.getNodeById(msc.getUserId()), I)); //run task on current node 
			double f1_func = Collections.max(f1); //f1 = MAX
			
			if(dag.outgoingEdgesOf(msc).isEmpty()) //if this task has no outgoing edges => finished with calculation 
				msc.setRank(f1_func);  //i.e. the current value is used as rank 
			else
			{ // if some additional dependent node exist (of the DAG)
								
				double tmpWRank;
				double maxSRank = 0; //maximum rank 
				for(ComponentLink neigh : dag.outgoingEdgesOf(msc)) // loop through all dependent tasks (via edge dependent) 
				{ // => i.e. all immediate succressors are investigated!
					tmpWRank = upRank(neigh.getTarget(),dag,I); // recursive: get rank of successor (needed for formula) 
					double f2_func = 0; // needed for computation of f2 value 
					List<Double> f2 = new ArrayList<Double>(); //same as for f1 
					if(neigh.getTarget().isOffloadable())  // if the successor can be offloaded (else the communication cost is 0!):
					{
						for(ComputationalNode cn : I.getAllNodes())  {//loop through nodes in I 
							f2.add(I.getTransmissionTime(neigh.getTarget(), I.getNodeById(msc.getUserId()), cn)); //get transmittion time (= communication cost) between those tasks 
						}
						f2_func = Collections.max(f2); //f2 = MAX
					}
					double tmpRank = tmpWRank + f2_func; //add both values (as needed for the formula) 
					maxSRank = (tmpRank > maxSRank)? tmpRank : maxSRank; //compute maximum  for the final rank 
					// if the current value exceeds the maximum => set this as new maximum
				}
				msc.setRank(f1_func + maxSRank); //f1 + f2 
			}
		}
		return msc.getRank(); //now every task received the rank that corresponds to its priority
	}
	
}
