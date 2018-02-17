package sa;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import cgramodel.CgraModel;
import cgramodel.PEModel;
import graph.CDFG;
import graph.Node;
import operator.Memory;
import operator.Operator;
import target.Processor;

public class Mutator {
	
	private static final int minNrPEs = 3;
	private static final int maxNrPEs = 16;
	private static final int nrDMAPEs = 4;
	
	private static final double interconnectRate = 0.015;
	private static final double operatorsRate    = 0.02;
	private static final double peCountRate	     = 0.2;
	
	
	private int nr;
	

	Set<Operator> necessaryOps;
	Set<Operator> DMAOps;
	
	public Mutator(CDFG graph) {
		necessaryOps = new HashSet<Operator>();
		DMAOps = new HashSet<Operator>();
		for (Node nd : graph.getNodes()) {
			necessaryOps.add(nd.getOperation());
			
		}
		for (Operator op : necessaryOps) {
			if (op.toString().contains("DMA")) {
				DMAOps.add(op);
			}
		}
		nr = 0;
	}// constructor
	
	public CgraModel applyMutation(CgraModel oldModel) {
		CgraModel newModel = null;
		
		System.out.println("mutator: mutate model #"+nr);
		
		applyMutationPECount(oldModel, peCountRate);
		applyMutationInterconnect(oldModel, interconnectRate);
		//applyMutationOperators(oldModel, operatorsRate);
		
		newModel = oldModel;
		newModel.finalizeCgra();
		
		nr++;
		return newModel;
	}// applyMutation
	
	
	private  void applyMutationPECount(CgraModel model, double rate) {
		double rnd = Math.random();
		if(rnd>=(1-rate)) {
			int currentCount = model.getNrOfPEs();
			double delta = Math.random()*4;
			delta -=2;
			int intDelta = (int)delta;
			int newNrPEs = model.getNrOfPEs()+intDelta;
			if (newNrPEs<3) {
				newNrPEs = 3;
			}
			if (newNrPEs>16) {
				newNrPEs = 16;
			}
			System.out.println("old Nr PEs = "+model.getNrOfPEs()+ "; new nr PEs = "+newNrPEs);
			//decrease if necessary
			while(newNrPEs<model.getNrOfPEs()) {
				int index = (int)(Math.random()*model.getNrOfPEs()) % model.getNrOfPEs();
				model.removePE(model.getPEs().get(index));
				
			}
			//increase if necessary
			while(newNrPEs<model.getNrOfPEs()) {
				model.addPE(getNewPE(model));
			}
			int minPEs = model.getNrOfPEs() < 4 ? 3 : 4;
			while (model.getNrOfMemoryAccessPEs()<minPEs ) {
				int index = (int)(Math.random()*model.getNrOfPEs()) % model.getNrOfPEs();
				PEModel pe = model.getPEs().get(index);
				for (Operator op : DMAOps) {
					pe.addOperator(op);
				}
				
			}
		}	//if rate fulfilled
	}//applyMutationPECount

	private  void applyMutationInterconnect(CgraModel model, double rate) {
		boolean[][] interconnect = new boolean[model.getNrOfPEs()][];
		for (int ii = 0; ii < model.getNrOfPEs(); ii++) {
			interconnect[ii] = new boolean[model.getNrOfPEs()];
			for (int jj = 0; jj < model.getNrOfPEs(); jj++) {
				if (model.getPEs().get(ii).getInputs().contains(model.getPEs().get(jj))) {
					interconnect[ii][jj] = true;
				}
			}
		}
		double rnd = Math.random();
		double interconnectThreshold = 0.988;
		if (rnd>(1-rate)) {
			for (int ii = 0; ii < model.getNrOfPEs(); ii++) {
				List<PEModel> inputsII = new LinkedList<PEModel>();
				for (int jj = 0; jj < model.getNrOfPEs(); jj++) {
					if (Math.random()>interconnectThreshold) {
						interconnect[ii][jj] = !(interconnect[ii][jj]);
						
					}//mutate input
					if (interconnect[ii][jj]) {
						inputsII.add(model.getPEs().get(jj));
					}// set input
				}//jj
				model.getPEs().get(ii).setInputs(inputsII);
			}//ii
		}//if random
		
	}//applyMutationInterconnect
	
	
	
	private  void applyMutationOperators(CgraModel model, double rate) {
		List<?> opList = Processor.Instance.getImplementedOperators();
		for (PEModel pe : model.getPEs()) {
			if (Math.random() > rate) { //mutate pe
				boolean[] opsInPE= new boolean[opList.size()];
				for (int ii = 0; ii < opsInPE.length; ii++) {
					if (pe.getAvailableOperators().containsKey(opList.get(ii))) {
						opsInPE[ii] = true;
						
					}
					if(Math.random() > 0.8) {
						opsInPE[ii] = !opsInPE[ii];
					}
					if(opsInPE[ii]) {
						pe.addOperator((Operator)opList.get(ii));
					} 
					else {
						if ((Operator)opList.get(ii)!=null) {
							pe.removeOperation((Operator)opList.get(ii));
						}
					}
				}// ii
			}//mutate
		}// pe
		int minPEs = model.getNrOfPEs() < 4 ? 3 : 4;
		while (model.getNrOfMemoryAccessPEs()<minPEs ) {
			int index = (int)(Math.random()*model.getNrOfPEs()) % model.getNrOfPEs();
			PEModel pe = model.getPEs().get(index);
			for (Operator op : DMAOps) {
				pe.addOperator(op);
			}
			System.out.println("stuck in while loop"+this.getClass());
		}
	}//applyMutationOperators
	
	private PEModel getNewPE(CgraModel model) {
		PEModel pe = new PEModel();
		int index = (int) Math.floor((Math.random()*123) % model.getNrOfPEs());
		PEModel peOld = model.getPEs().get(index);
		for (Operator op : peOld.getAvailableOperators().keySet()) {
			pe.addOperator(op, peOld.getAvailableOperators().get(op));
		}
		for (PEModel pe2 : model.getPEs()) {
			pe.addPE2inputs(pe2);
		}
		pe.setLiveout(true);
		pe.setLivein(true);
		return pe;
	}//


}
