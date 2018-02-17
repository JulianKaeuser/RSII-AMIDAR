package sa;

/**
 * 
 */


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import graph.CDFG;

import graph.LG;
import graph.Loop;
import graph.Node;

import scheduler.Interval;
import scheduler.MissingOperationException;
import scheduler.NotSchedulableException;
import scheduler.Schedule;
import scheduler.Scheduler;
import target.Amidar;
import target.Processor;


/**
 * @author juliank
 *	This class is s <Scheduler> which does an ASAP scheduling, involving nested loops and doing node-fusing
 *	where possible. The schedule is independent of any hardware and only displays the ASAP version of
 *	the graph. This version may be used to estimate e.g. the minimum schedule length any other, resource- or
 *	time constrained <Scheduler> can achieve. It may also serve as a initial step for other schedulers
 */
public class ASAPScheduler implements Scheduler {
	public static final boolean debug = true;
	public static final boolean consolePrint = true;
	public static final  boolean dumpSchedules = true;
	
	/** The graph to estimate the schedule for */
	private CDFG cdfg;
	private LG lg;
	
	
	/** The schedule produced */
	private Schedule schedule;
	
	
	/**
	 * Constructor which initially sets the graph
	 * @param rg
	 */
	public ASAPScheduler(){

	}
	
	


	
	/**
	 * @return the schedule
	 */
	public Schedule getSchedule() {
		return schedule;
	}



	
	/**
	 * Determines the lower bound for the schedule length (unconstrained ASAP schedule length)
	 * for the graph stored in this ScheduleEstimation instance
	 * @return the length of the unconstrained ASAP schedule; if errors occur, e.g. the graph==null,
	 *  -1 is returned
	 */
	public int getMinimumScheduleLength() throws NotSchedulableException{
		if (cdfg ==null) {
			throw new NotSchedulableException();
		}
		schedule = new Schedule();
		

		
		/** holds all scheduled nodes*/
		Set<Node> scheduled = new HashSet<Node>();
		
		/** holds all candidates for the current time step*/
		Set<Node> candidates = new HashSet<Node>();
		
		/** holds all unscheduled nodes*/
		Set<Node> unscheduled = new HashSet<Node>(cdfg.getNodes());
		
		/** holds all nodes which have been scheduled in the (at least) previous timestep and might be unfinished*/
		Set<Node> pending = new HashSet<Node>();
		
		/** Holds all nodes which must be removed from the candidate set*/
		Set<Node> rmFromCandidates = new HashSet<Node>();
		
		/** Maps the starting times and ending time steps (context-related) for each loop
		 * 	[0] = start
		 *  [1] = end
		 * */
		Map<Loop, int[]> loopTimes = new HashMap<Loop, int[]>();
		
		/** The loop which is currently being scheduled*/
		Loop currentLp = lg.getRoot();
		if (currentLp == null) {
			System.err.println("Loopgraph has no root!");
			throw new NotSchedulableException();
		}
		
		/** the current time step which is planned*/
		int currentTs = 0;
		
		Set<Node> electorates = new HashSet<Node>();
		for (Node nd : unscheduled){
			if (cdfg.getAllPredecessors(nd).isEmpty()){
				if (currentLp.contains(nd)){
					// node can be scheduled
					scheduled.add(nd);
					electorates.add(nd);
					if (isNonNative(nd)){
						schedule.add(nd, new Interval(currentTs, currentTs+getLatency(nd)));
					}
				}
			}
		}
		unscheduled.removeAll(electorates);
		electorates.clear();
		
		// mark all constants scheduled
		for (Node nd : cdfg.getNodes()){
			if (nd.getOperation().isConst()){
				scheduled.add(nd);
				unscheduled.remove(nd);
			}
		}
		
		printCurrentSchedule();
		
		while (!unscheduled.isEmpty()){
			candidates.clear();
			rmFromCandidates.clear();
			Set<Node> rmFromPending = new HashSet<Node>();
			for (Node nd : pending){
				if (schedule.slot(nd).ubound<=currentTs){
					rmFromPending.add(nd);
				}
				for (Node storeNd : cdfg.getAllSuccessors(nd)){
					if (storeNd.getOperation().isRegfileStore()){
						scheduled.add(storeNd);
						unscheduled.remove(storeNd);
					}
				}
			}
			pending.removeAll(rmFromPending);
			for (Node nd : unscheduled){
				/*
				 * General compatibility for Node nd is checked:
				 * is candidate if:
				 * 1 - all predecessors are scheduled
				 * 2 - all predecessors are finished
				 */
				if (isCandidate(nd, cdfg, scheduled, pending)){
					candidates.add(nd);
				}
			}//end nd: unscheduled
			//check for fuseable load nodes
			for (Node nd : candidates){
				if (nd.getOperation().isRegfileLoad()){
					//fuse consumer of load node and 
					scheduled.add(nd);
					unscheduled.remove(nd);
					for (Node ndSucc : cdfg.getAllSuccessors(nd)){
						if (isCandidate(ndSucc, cdfg, scheduled, pending)){
							candidates.add(ndSucc);
						}
					}
					
				}
			}// fusing of load nodes done
			//place candidates in scheduled, in schedule, in pending
			Set<Node> rmCands = new HashSet<Node>();
			for (Node nd : candidates){
				Loop ndLoop = lg.getLoop(nd);
				Loop pivotLp = checkLoop(nd, scheduled, ndLoop, currentLp, currentTs, pending); 
				if (pivotLp!=null){
					currentLp = pivotLp;
					scheduled.add(nd);
					unscheduled.remove(nd);
					rmCands.add(nd);
					// add only real nodes to schedule and pending
					if (isNonNative(nd)){	
						Interval slot = new Interval(currentTs, currentTs+getLatency(nd)-1);
						schedule.add(nd, slot);
						pending.add(nd);
					}
				
				}//endif checkLoop
			}
			candidates.removeAll(rmCands);
			currentTs++;
		}
		
		return schedule.length();
	}
	
	/**
	 * Returns the latency connected to this Node's Operator
	 * @param nd the node to get the latency from
	 * @return the latency; if the Operator is native, 0 is returned (should not be scheduled anyway)
	 */
	private int getLatency (Node nd){
		return nd.getOperation().getImplementation().getLatency();
	}
	
	/**
	 * Reveals whether this node can be scheduled and is not a native operation
	 * (STORE, LOAD, CONST, CONST64, STORE64, LOAD64, NOP, ERROR cannot be scheduled)
	 * @param nd the Node to check
	 * @return true if the operator is not from the list, false if it is
	 */
	private boolean isNonNative(Node nd){
		if (nd.getOperation().isNative()) {
			return true;
		}
		return true;
	}
	
	/**
	 * Prints the current schedule
	 * @return a String with the schedule
	 */
	private String printCurrentSchedule (){
		String str = schedule.toString();
		if (consolePrint){
			System.out.println(str);
		}
		return str;
	}
	

	
	/**
	 * Checks whether the given node may be scheduled, depending on
	 * its predecessors and with regard to the loops
	 * @param nd the Node to check
	 * @param scheduled the Set of scheduled nodes
	 * @param ndLoop the loop of the node nd
	 * @return 
	 * null if not possible to schedule the node
	 * the loop in which the scheduler has to go on (new currentLoop) if the node can be scheduled
	 */
	private Loop checkLoop(Node nd, Set<Node> scheduled, Loop ndLoop, Loop currentLoop, int currentTs, Set<Node> pending){
		
		Set<Node> nodesOfCurrentLoop = new HashSet<Node>();
		Set<Node> nodesOfNdLoop = new HashSet<Node>();
		for (Node node : cdfg.getNodes()){
			if (lg.getLoop(node).equals(currentLoop)){
				nodesOfCurrentLoop.add(node);
			}
		}
		for (Node node : cdfg.getNodes()){
			if (lg.getLoop(node).equals(ndLoop)){
				nodesOfNdLoop.add(node);
			}
		}
		
		if (ndLoop.equals(currentLoop)){
			//if no other nodes from other loops than ndLoop have been scheduled in this ts
			for (Node ndOther : schedule.nodes(currentTs)){
				if (!lg.getLoop(ndOther).equals(ndLoop)){
					return null; //cant schedule node: at least one node scheduled in this ts is from another loop
				}
			}
			return currentLoop; //can schedule
		}//end currentLoop-> currentLoop transition
		
		if (lg.isChildOf(ndLoop, currentLoop)){
			//candidate nodes' loop is child of the current
			// can only be scheduled if
			//(1) no nodes from superloop (currentLoop) are scheduled in the current ts
			//(2) no nodes are pending
			//(3) all predecessors necessary for the subloop(ndLoop) are finished(scheduled& not pending)
			for (Node ndOther : schedule.nodes(currentTs)){
				if (!lg.getLoop(ndOther).equals(ndOther)){
					return null; //case (1): at least one node scheduled in this ts is drom antoher loop
				}
			}
			if (!pending.isEmpty()){
				return null; //case (2): pending nodes; can't begin new subloop
			}
			for (Node ndOther : nodesOfNdLoop){	//for each node in the new subloop
				for (Node ndOtherPred : cdfg.getAllPredecessors(ndOther)){ //for each predecessor
					if (!lg.getLoop(ndOtherPred).equals(ndLoop)){ 	//check if any predecessor that is not from new subloop
						if (!scheduled.contains(ndOtherPred) || pending.contains(ndOtherPred)){ //is not finished
							return null; //case (4): at least one result necessary for one node in the new subloop (ndLoop) but coming from another lop is not finished
						}
					}
				}
			}
			return ndLoop; //can schedule
		}//end superloop->subloop transition
		
		if (lg.isChildOf(currentLoop, ndLoop)){
			//candidate nodes' loop is parent of the current
			//can only be scheduled if 
			// - no nodes from the currentLoop have been scheduled in this ts (1)
			// -no nodes from the currentLoop are pending (2)
			// all nodes of currentLoop are scheduled
			for (Node ndOther : schedule.nodes(currentTs)){
				if (lg.getLoop(ndOther).equals(currentLoop)){
					return null; //case (1): cant be scheduled(=go back to superloop): 
								// there is at least one node scheduled in the currentTs which has the currentLoop and not another
				}
			}
			for (Node ndOther : pending){
				if (currentLoop.contains(ndOther)){
					return null; //case (2): cant be scheduled(=go back to superLoop) because there are pending nodes from the currentLoop
				}
			}
			for (Node ndOther : nodesOfCurrentLoop){
				if (!scheduled.contains(ndOther)){
					return null; //case (3): at least one node from the currentLoop (subloop) is not scheduled
				}
			}
			return ndLoop;// can schedule
		}// end subloop->superloop transition
		
		if (lg.getFather(ndLoop).equals(lg.getFather(currentLoop))&& !currentLoop.equals(ndLoop)) { //both are subloops on the same level
			//can only be scheduled if
			//(1) no nodes from other loop are scheduled in the current ts
			//(2) no nodes of the last subloop are pending
			//(3) all nodes of the last subloop are scheduled
			//(4) all pred. nodes necessary for the next subloop (ndLoop) are finished
			for (Node ndOther : schedule.nodes(currentTs)){
				if (!lg.getLoop(ndOther).equals(ndLoop)){
					return null; //case (1): at least one other node from another loop than the one planned is yet scheduled
				}
			}
			for (Node ndOther : pending){
				if (currentLoop.contains(ndOther)){
					return null;//case (2): at least one other node from the previous same-level loop is pending
				}
			}
			for (Node ndOther : nodesOfCurrentLoop){
				if (!scheduled.contains(ndOther)){
					return null; //case (3): at least one node from the currentLoop (subloop) is not scheduled
				}
			}
			for (Node ndOther : nodesOfNdLoop){	//for each node in the new subloop
				for (Node ndOtherPred : cdfg.getAllPredecessors(ndOther)){ //for each predecessor
					if (!lg.getLoop(ndOtherPred).equals(ndLoop)){ 	//check if any predecessor is not from new subloop
						if (!scheduled.contains(ndOtherPred) || pending.contains(ndOtherPred)){
							return null; //case (4): at least one result necessary for one node in the new subloop (ndLoop) but coming from another lop is not finished
						}
					}
				}
			}
			return ndLoop;
		}// end subloop-subloop (same level) transition
		return null;
	}
	
	/**
	 * Reveals whether the given Node nd is a schedule candidate, depending on the yet scheduled and pending
	 * nodes
	 * @param nd the Node to check
	 * @param dcfg the graph
	 * @param scheduled a Set of yet scheduled nodes
	 * @param pending s Set of pending nodes
	 * @return true if the node is a candidate, false if not
	 */
	private boolean isCandidate(Node nd, CDFG dcfg, Set<Node> scheduled, Set<Node> pending){
		Set<Node> ndPreds = dcfg.getAllPredecessors(nd);
		boolean isCandidate = true;
		for (Node ndPred : ndPreds){
			if (!scheduled.contains(ndPred)){
				isCandidate = false;
			} 
			if (pending.contains(ndPred)){
				isCandidate = false;
			}	
		}
		return isCandidate;
	}



	public Schedule schedule() throws NotSchedulableException,
			MissingOperationException {
		// TODO Auto-generated method stub
		if (cdfg== null || lg == null){
			throw new NotSchedulableException();
		}
		getMinimumScheduleLength();
		
		return schedule;
	}
	
	public String drawOrderSchedule(){
		String res = null;
		StringBuilder str = new StringBuilder();
		
		String dNSize = "!\", height = \"1.5\", width=\"1.5\"];\n ";
		String dNShape = "[shape=\"circle\", style=\"filled\", color=\"004E8ABF\", pos=\"";
		int maxNrNodes = 0;
		for (int ii = 0; ii<schedule.length(); ii++){
			int slotNodes = schedule.nodes(ii).size();
			if (slotNodes >=maxNrNodes){
				maxNrNodes = slotNodes;
			}
		}
		HashMap<Node, Integer> nodeWidth = new HashMap<Node, Integer>();	
		HashMap<Integer, Node[]> usedSlotWidths = new HashMap<Integer, Node[]>();
		for (int aa = 0; aa<schedule.length(); aa++){
			usedSlotWidths.put(aa, new Node[maxNrNodes]);
		}
		for (int ii= 0; ii<schedule.length(); ii++){
			for (Node nd : schedule.nodes(ii)){
				if (nodeWidth.get(nd)==null){
					int value = -1;
					for (int jj = 0; jj<usedSlotWidths.get(ii).length; jj++){
						if (usedSlotWidths.get(ii)[jj]==null){
							usedSlotWidths.get(ii)[jj] = nd;
							for (int bb = 1; bb<getLatency(nd); bb++){
								if (usedSlotWidths.get(ii+bb)!= null){
									usedSlotWidths.get(ii+bb)[jj] = nd;
								}
							}
							value = jj;
							break;
						}
					}
					nodeWidth.put(nd, value);	
				}
			}
		}
		
		int height = 0;
		int width = 0;
		double maxHeight = schedule.length()*2.5;
		str.append("//do not use DOT to generate pdf use NEATO or FDP\n ");
		str.append("digraph{\n ");
		str.append("splines=\"ortho\";\n ");
		
		while (height!=schedule.length()){
			height++;
			double currentPosHeight = maxHeight-(height*2.5);
			for (Node nd : schedule.nodes(height)){
				String ndName = "\""+nd.toString()+"\"";
				str.append(ndName);
				str.append(dNShape);
				Double ndHeight = currentPosHeight;
				Double ndWidth =  2.5*nodeWidth.get(nd);
				String pos = ndWidth.toString()+","+ndHeight.toString();
				str.append(pos);
				str.append(dNSize);
				Set<Node> ndSuccs = getFusedNdSuccs(nd);
				if (!(ndSuccs.isEmpty() || ndSuccs==null)){
					for (Node ndSucc : ndSuccs){
						str.append(ndName);
						str.append(" -> ");
						String ndSuccName = "\""+ndSucc.toString()+"\"";
						str.append(ndSuccName);
						str.append(";\n ");
					}
					str.append("\n ");
				}
			}
			
		}
		res = str.toString();
		return res;
	}

	/**
	 * Returns all direct successors (data dependencies) of the given node
	 * @param nd
	 * @return
	 */
	private Set<Node> getFusedNdSuccs(Node nd) {
		Set<Node> ndSuccs = new HashSet<Node>();
		if (cdfg.getConsumers(nd)==null){
			return ndSuccs;
		}
		else{
		for (Node ndOther : cdfg.getConsumers(nd)){
			if (isNonNative(ndOther)){
				ndSuccs.add(ndOther);
			}
			if (ndOther.getOperation().isRegfileAccess()){
				if (cdfg.getConsumers(ndOther)!=null){
					ndSuccs.addAll(cdfg.getConsumers(ndOther));
				}
			}
		}
		}
		return ndSuccs;
	}
	
	private void printSchedule(String path){
		File file = new File(path);
		FileWriter writer;
		try {
			writer = new FileWriter(file);
			BufferedWriter buf = new BufferedWriter(writer);
			for (int ii= 0; ii<schedule.length(); ii++){
				buf.write(ii+";");
				for (Node nd : schedule.nodes(ii)){
					buf.write(nd.toString()+";");
				}
				buf.newLine();
				buf.flush();
			}
			buf.flush();
			buf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}


	public void setGraphs(CDFG graph, LG lg) {
		// TODO Auto-generated method stub
		this.cdfg = graph;
		this.lg = lg;
	}
}