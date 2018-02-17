package sa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import cgramodel.CgraModel;
import graph.CDFG;
import io.AttributeParserAmidar;
import io.AttributeWriterAmidar;
import scheduler.Schedule;

public class SimulatedAnnealing {
	
	private static final double innerNum = 7;
	
	private static final double nPow43 = 2100; //approx. 16 (PEs) * 20 (ops) ^4/3
	
	private String pathToCurrentOptimumConfig;
	
	private CgraModel bestModel;
	
	private Schedule schedule;
	
	private CostFunction cf;
	
	private Mutator mutator;
	
	private String configDir;
	
	private CgraModel currentConfig;
	
	//io
	private AttributeParserAmidar aParser;
	private AttributeWriterAmidar aWriter;
	
	//sa algorithm stuff
	private int innerRuns;
	private double currentTemperature;
	private double acceptanceRate;
	
	private static double startingCost = Double.MAX_VALUE;
	
	public SimulatedAnnealing(CostFunction costfunction, String cgraConfigDir, CgraModel startingComposition, CDFG cdfg) {
		cf = costfunction;
		configDir = cgraConfigDir;
		currentConfig = startingComposition;
		bestModel = startingComposition;
		aParser = new AttributeParserAmidar();
		aWriter = new AttributeWriterAmidar();
		
		mutator = new Mutator(cdfg);
		
		innerRuns = 0;
		acceptanceRate = 1;
		
	}//constructor
	
	public void run() {
		// perform the simulated annealing main loop here
		/*
		S=RandomConfiguration();
		T=InitialTemperature();
		while (ExitCriterion()==false) {
			while (InnerLoopCriterion() == false) {
				S new = Generate(S);
				ΔC = Cost(S new )-Cost(S);
				r = random(0,1);
				if (r < e -ΔC/T ){
				 	S = S new
				}
			}
			T=UpdateTemperature();
		}
		*/
		currentTemperature = initTemp();
		double currentCost = cf.getCost(currentConfig);
		startingCost = currentCost;
		CgraModel oldModel = currentConfig;
		
		double currentBestCost = Double.MAX_VALUE;
		
		int iterationCounter = 1;
		int accepted = 0;
		while (!outerCriterion()) {
			System.out.println("SA current Temperature: "+currentTemperature);
			while(!innerCriterion()) {
				String tempDir = configDir+ "temp/";
				aWriter.writeJSON(oldModel, "cgra_temp", tempDir);
				CgraModel temp = aParser.loadCgra(tempDir+"cgra_temp.json");
				CgraModel newModel = mutator.applyMutation(temp);
				double newCost = cf.getCost(newModel);
				System.out.println("nc:"+newCost);
				if (newCost < currentBestCost) {
					currentBestCost = newCost;
					aWriter.writeJSON(newModel, "best", configDir+"best/");
					System.out.println("sa best cost: "+currentBestCost/(10240));
					dump(newModel, currentBestCost,iterationCounter, configDir+"best/");
				}
				//System.out.println("cost of current sol: "+newCost);
				double deltaCost = newCost-currentCost;
				
				double r = cf.isValid() ? Math.random() : 1.01;
				double quot = -deltaCost / currentTemperature;
				double ePowerQuot = Math.pow(Math.E, quot);
				if (r < ePowerQuot) {
					oldModel = newModel;
					currentCost = newCost;
					accepted ++;
					acceptanceRate = accepted/iterationCounter;
				}// accept
				iterationCounter++;
			}//inner
		}//outer
		
	}//
	
	private boolean outerCriterion() {
		if (currentTemperature< 10) {
			currentTemperature = initTemp();
			return true;
		}
		double gamma = 0.8;
		if (acceptanceRate>0.15) {
			gamma = 0.95;
		}
		if (acceptanceRate > 0.8) {
			gamma = 0.9;
			if (acceptanceRate>0.96) {
				gamma = 0.5;
			}
		}
		
		currentTemperature *= gamma;
		return false;
	}// outerCriterion
	
	private boolean innerCriterion() {
		innerRuns++;
		if (innerRuns > nPow43) {
			innerRuns=0;
			return true;
		}
		return false;
	}//innerCriterion
	
	public double initTemp() {
		return 8000.0;
	}
	
	public static String returnResultConfigurationPath() {
		
		return null;
	}
	
	public static void dump(CgraModel model, double cost, int iteration, String path) {
		try {
			BufferedWriter buf = new BufferedWriter(new FileWriter(new File(path+"info.txt")));
			buf.append("cost: "+cost);
			buf.newLine();
			buf.append("starting cost was: "+startingCost);
			buf.append("iteration nr. "+iteration);
			buf.newLine();
			buf.append("############################################");
			buf.newLine();
			buf.append(model.toString());
			buf.newLine();
			buf.flush();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
