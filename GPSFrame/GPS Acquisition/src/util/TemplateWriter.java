package util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.stringtemplate.v4.ST;
import org.stringtemplate.v4.STGroup;
import org.stringtemplate.v4.STGroupFile;

public class TemplateWriter {
	
	private final static int samplingRate = 400000; //400kHz
	private final static float halfSamplingRate = (float)0.5 * (float)samplingRate;
	private final static float oneDivSamplingRate = 1/samplingRate;
	private final static int stepRate 	= 1000; //1 kHz
	private final static int maxRate	 	= 5000; //5kHz
	private final static int minRate 		= -5000; //-5kHz
	
	private final static int nrFrequencies = Math.floorDiv((maxRate-minRate), stepRate)+1;
	
	private final static int[] allFrequencies ={-5000, -4000, -3000, -2000, -1000, 0, 1000, 2000, 3000, 4000, 5000};
	
	private final static float[] allFrequenciesDivSamplingRate = null;
	
	private final static float[][][] precomputedFrequencyShifts = null;
	
	private final static float gamma 		= (float)0.015;
	
	private final static int nDivisor		= 16;
	
	private final static int maxFactorPreCompute = 256;
	
	private  int nrOfSamples;
	private  float oneDivN;

	public static void main (String[] args) {
		String javaFileContent = null;
		String outFileNameTemplate = "AqcuisitionTemplate.stg";
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File("src/Aqcuisition.java")));
			BufferedWriter buf = new BufferedWriter(new FileWriter(new File("AqcuisitionTempalte.stg")));
			String line;
			buf.write("Aqcuisition ( preComputedFrequencies\n");
			buf.write(")\n ::=<<\n");
			while ((line = reader.readLine())!=null) {
				if (line.startsWith("private final static float[][] precomputedFrequencyShiftsReal =")) {
					line += "§precomputedFrequencyShifts§";
				}
				
				buf.write(line);
			}
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		//2nd step: write into newly generated Template
		STGroupFile stgroup = new STGroupFile("AqcuisitionTemplate.stg", '§', '§');
		ST template = stgroup.getInstanceOf("Acquisition");
		
		
		try {
			BufferedWriter buf = new BufferedWriter(new FileWriter(new File("Acquisition.java")));
			buf.write(template.render());
			buf.flush();
			buf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static String printAllFrequenciesDivSamplingRateToConsole() {
		StringBuilder b = new StringBuilder();
				
				b.append("{");
					int[] f = new int[nrFrequencies];
					for (int fd = 0; fd < nrFrequencies; fd++) {
						if (fd!=0) {
							b.append(", \n");
						}
						float val = ((float)minRate+((float)fd*(float)stepRate))*(float)Math.PI
								/(float)halfSamplingRate;
						b.append("(float)"+new Float( val).toString());	
					}	
				b.append("\n};");
				return b.toString();
			}
}
