package sa;

import cgramodel.CgraModel;
import cgramodel.PEModel;
import graph.CDFG;
import graph.LG;
import graph.Loop;
import operator.Implementation;
import scheduler.MissingOperationException;
import scheduler.NotSchedulableException;
import scheduler.RCListSched;
import scheduler.Schedule;

public class CostFunctionEnergy implements CostFunction{

	private static final int N = 400;
	private static final int stages = 10;
	
	private boolean isValid;
	
	private CDFG cdfg;
	private LG lg;
	private int localVars;
	
	public CostFunctionEnergy (CDFG graph, LG loopGraph, int maxNumLocalVars) {
		localVars = maxNumLocalVars;
		cdfg = graph;
		lg = loopGraph;
	}//constructor

	/**
	 * Returns the cost for this cgramodel and the given graph/lg combination.
	 * Changes the isValid field if necessary
	 * @param model
	 * @param cdfg
	 * @param lg
	 * @return
	 */
	public double getCost(CgraModel model) {
		double cost = Double.MAX_VALUE;
		//TODO implement cost function for simulated annealing
		
		try {
			RCListSched scheduler = new RCListSched(cdfg, lg, localVars, model);
			Schedule schedule = scheduler.schedule();
			int length = schedule.length();
			System.out.println("schedlength = "+length);
			
			int[] loopLengths = new int[lg.getMaxNestingLevel()+1];
			int[] nestingLevels = new int[lg.getMaxNestingLevel()+1];
			int counter = 0;
			for (Loop loop : lg.loops) {
				int[] times = scheduler.getLoopTimes().get(loop);
				loopLengths[counter] = times[1] - times[0]; 
				
				nestingLevels[counter] = lg.getNestingLevel(loop);
				System.out.println("cf: loop "+loop+"; nesti level: "+nestingLevels[counter]+ "; len: "+(times[1]-times[0]));
				counter++;
			}//loop
			int maxNestingLevel = 0;
			int maxIndex = 0;
			int maxLoopLength = 0;
			for (int ii = 0; ii < loopLengths.length; ii++) {
				if (nestingLevels[ii]>= maxNestingLevel) {
					maxNestingLevel = nestingLevels[ii];
					
					maxIndex = ii;
					maxLoopLength = loopLengths[ii];
					//System.out.println("cost function: nesting level "+ maxNestingLevel + ", length "+maxLoopLength);
					
				}
			}
			cost = 1024*10*maxLoopLength;
			
			// energy consumption
			cost = cost*getEnergy(model);
			isValid = true;
			
			
		} catch (MissingOperationException e) {
			// TODO Auto-generated catch block
			System.out.println("not scheduleable: missing op");
			e.printStackTrace();
			
			isValid = false;
			return cost;
		} catch (NotSchedulableException e) {
			// TODO Auto-generated catch block
			System.out.println("not scheduleable: ");
			e.printStackTrace();
			isValid = false;
			return cost;
		}
		
		return cost;
	}

	public boolean isValid() {
		// TODO Auto-generated method stub
		return isValid;
	}
	
	private static double getEnergy(CgraModel model) {
		double en = 1.0;
		
		for (PEModel pe : model.getPEs()) {
			for (Implementation impl : pe.getAvailableOperatorImplementations()) {
				en += impl.getEnergyconsumption();
			}
		}
		return en;
	}
}
