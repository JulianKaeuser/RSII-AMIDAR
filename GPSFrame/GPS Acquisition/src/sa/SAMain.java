package sa;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import cgramodel.CgraModel;
import graph.CDFG;
import graph.LG;
import io.AttributeParserAmidar;
import io.AttributeWriterAmidar;
import io.ObjectSerializer;
import scheduler.NotSchedulableException;

public class SAMain {

	public static void main(String[] args) {
		String cgraConfigDir = args[0];
		String startingConfig = args[1];
		String cdfgPath = args[2];
		String lgPath = args[3];
		String localVarsPath = args[4];
		
		System.out.println("target dir: "+cgraConfigDir);
		System.out.println("startingConfig: "+startingConfig);
		System.out.println("cdfg: \t   "+cdfgPath);
		System.out.println("lg: \t   "+lgPath);
		System.out.println("localVars: "+localVarsPath);
		
		ObjectSerializer deser = new ObjectSerializer();
		int numLocalVars = findNumLocalVars(localVarsPath);
		System.out.println("localVars value: "+numLocalVars);
		
		
		//LG loopgraph = (LG) deser.deserialize(lgPath);
		//System.out.println("LG Value: "+loopgraph.toString());
		//CDFG graph = (CDFG) deser.deserialize(cdfgPath);
		
		String[] leftStrings = {args[0], args[1]};
		CDFG graph = null;
		LG loopGraph = null;
		int maxNumLocalVars = numLocalVars;
		runSA(graph, loopGraph, maxNumLocalVars, leftStrings);
		

	}//main
	
	public static void runSA(CDFG graph, LG loopGraph, int maxNumLocalVars, String[] args) {
		String cgraConfigDir;
		String startingConfig;
		if (args==null) {
			cgraConfigDir = "/home/juliankaeuser/Dokumente/Studium2/Rechnersysteme2/Uebungen/Uebung2/AmidarRS2/amidar-sim2/GPSFrame/GPS Acquisition/sa_target/";
			startingConfig = "/home/juliankaeuser/Dokumente/Studium2/Rechnersysteme2/Uebungen/Uebung2/AmidarRS2/amidar-sim2/Amidar/config/FU/CGRA/CGRA_16.json";
		}
		else {
		 cgraConfigDir = args[0];
		 startingConfig = args[1];
		}
		
		
		System.out.println("target dir: "+cgraConfigDir);
		System.out.println("startingConfig: "+startingConfig);
		
		
		
		//cgra
		AttributeParserAmidar modelParser = new AttributeParserAmidar();
		CgraModel startingComposition = modelParser.loadCgra(startingConfig);
		startingComposition.setName("test");
		
		AttributeWriterAmidar modelWriter = new AttributeWriterAmidar();
		modelWriter.writeJSON(startingComposition, "test", cgraConfigDir);
		
		CostFunction cf = new CostFunctionEnergy(graph, loopGraph, maxNumLocalVars);
		SimulatedAnnealing sa = new SimulatedAnnealing(cf, cgraConfigDir, startingComposition, graph);
		
		/*ASAPScheduler asap = new ASAPScheduler();
		//asap.setGraphs(graph, loopGraph);
		
		try {
			int minLength = asap.getMinimumScheduleLength();
			System.out.println("min length: "+minLength);
		} catch (NotSchedulableException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
		sa.run();
		System.out.println("finished SA");
	}
	
	private static int findNumLocalVars(String path) {
		JSONParser parser;
		JSONObject json = null;
		parser = new JSONParser();
		FileReader fileReader;
		
		try {
			fileReader = new FileReader(path);
			json = (JSONObject) parser.parse(fileReader);
			
			long maxValue = (Long) json.get("value");
			return new Integer( new Long(maxValue).intValue()).intValue();
		} 
		catch (FileNotFoundException e) {
			System.err.println("No file found- SAMain");
			e.printStackTrace(System.err);
	
		} 
		catch (IOException e) {
			System.err.println("Error while reading file ");
			e.printStackTrace(System.err);
		} 
		catch (ParseException e) {
			System.err.println("Error while reading file ");
			e.printStackTrace(System.err);
		}
		return 0;
	}//finNumLocalVars

}
