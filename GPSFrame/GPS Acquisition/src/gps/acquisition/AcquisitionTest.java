package gps.acquisition;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;


public class AcquisitionTest {

	public static void main(String[] args) {
		

		AcquisitionNoOpt acq = null;
		
		
		int testNr = 0;
		while(testNr <1 || testNr >4) {

			 testNr =  (int)Math.round(Math.random()*(4.0));
		}
		
		String Nr = Integer.toString(testNr);
		if(Nr.length() == 1){
			Nr ="0"+ Nr;
		}
		
		String inputData = "../Amidar/input/"+Nr+"_inputdata.dat";
		String inputCodes = "../Amidar/input/"+Nr+"_inputcodes.dat";
		String results = "../Amidar/input/"+Nr+"_results.dat";
		
		boolean acquisition = false;
		int freq = -1234;
		int codeVersch = -1110;
		int nrOfSamples = -1;

		try{
			FileInputStream frData = new FileInputStream(inputData);
			BufferedInputStream brData = new BufferedInputStream(frData);
			
			
			
			String line = readLine(brData);
			nrOfSamples = Integer.parseInt(line);
			acq = new AcquisitionNoOpt(nrOfSamples);
			
			
			for(int i = 0; i < nrOfSamples; i++){
				line = readLine(brData);
				line = line.substring(1, line.length()-1);
				String [] values = line.split(",");
				float real = Float.parseFloat(values[0]);
				float imag = Float.parseFloat(values[1]);
				acq.enterSample(real, imag);
				
			}

			brData.close();
			
			FileInputStream frCode = new FileInputStream(inputCodes);
			BufferedInputStream brCode = new BufferedInputStream(frCode);
			
			line = readLine(brCode);
			if(nrOfSamples != Integer.parseInt(line)){
				System.out.println("ERROR: Code length and data length do not match");
				System.exit(0);
			}
			

			for(int i = 0; i < nrOfSamples; i++){
				line = readLine(brCode);
				line = line.substring(1, line.length()-1);
				String [] values = line.split(",");
				float real = Float.parseFloat(values[0]);
				float imag = Float.parseFloat(values[1]);
				acq.enterCode(real, imag);
				
			}
			
			FileInputStream fsRes = new FileInputStream(results);
			BufferedInputStream bsRes = new BufferedInputStream(fsRes);
			
			

			line = readLine(bsRes);
			
			acquisition = line.split(" ")[0].equals("positive");
			
			codeVersch = Integer.parseInt(readLine(bsRes).split("=")[1].trim());
			freq = Integer.parseInt(readLine(bsRes).split("=")[1].trim());


		} catch(IOException e){
			System.err.println("Error while reading input data");
			e.printStackTrace(System.err);

		}

		
		boolean res = acq.startAcquisition();
		
		
			boolean passed = res == acquisition && acq.getCodeVerschiebung()== codeVersch && acq.getDopplerverschiebung() == freq;
			

		System.out.println("[Test Results]");
		System.out.println((passed?"PASSED":"FAILED")+ " Test Nr. " + Nr);
		
			System.out.println("nrOfSamples: "+nrOfSamples);
			System.out.println("Expected " + acquisition + " acquistion, got "+acquisition);
			System.out.println("    codeShift:    " + codeVersch + ", got: "+acq.getCodeVerschiebung());
			System.out.println("    dopplerShift: " + freq+ ", got: "+acq.getDopplerverschiebung());
		
		


	}
	
	private static String readLine(BufferedInputStream bs){
		char[] ret = new char[200];
		int cnt = 0;
		try{
		int val = bs.read();
		
		while(val != '\n'){
			ret[cnt++] = (char)val;
			val = bs.read();
		}
		} catch (IOException e){
			
		}
		
		return new String(ret, 0, cnt);
		
		
	}

}
