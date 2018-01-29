package util;

public class Util {
	
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
		//printPrecomputedFrequenciesToConsole();
		printAllFrequenciesDivSamplingRateToConsole();
		//printRateStepToConsole();
		//printPrecomputedFrequencyShiftsToConsole(true);
	}
	
	public static void printRateStepToConsole() {
		int rate = Math.floorDiv((maxRate-minRate), stepRate)+1;
		System.out.println(rate);
	}
	
	public static void printPrecomputedFrequenciesToConsole() {
		StringBuilder b = new StringBuilder();
		
		b.append("{");
			int[] f = new int[nrFrequencies];
			for (int fd = 0; fd < nrFrequencies; fd++) {
				if (fd!=0) {
					b.append(", ");
				}
				b.append(new Integer( minRate+(fd*stepRate)).toString());	
			}	
		b.append("}");
		System.out.println(b.toString());
	}


	public static void printAllFrequenciesDivSamplingRateToConsole() {
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
		System.out.println(b.toString());
	}
	
	public static void printPrecomputedFrequencyShiftsToConsole(boolean real) {
		int index = real? 0 : 1;

			StringBuilder b = new StringBuilder();
				if(real) {
					b.append("real:{");
				}
				else {
					b.append("imag:{"); 
				}
				
				float[][][] shifts = new float[nrFrequencies][][];
				for (int fd = 0; fd< nrFrequencies; fd++) {
					if(fd!=0) {
						b.append(",\n");
						}
					b.append("{");
					shifts[fd] = new float[nDivisor*maxFactorPreCompute][];
					for (int nn = 0; nn < nDivisor*maxFactorPreCompute; nn++) {
						if(nn!=0) b.append(", \n");
						shifts[fd][nn] = getCartesianFormFromPolar(1, (float)-2*(float)Math.PI*nn*allFrequencies[fd]*oneDivSamplingRate);
						b.append(shifts[fd][nn][index]);
					} // n
					b.append("}");
				} // m
		
				b.append("\n}");
			System.out.println(b.toString());
		}
	
	public static float[] getCartesianFormFromPolar(float abs, float phase) {
		float[] res = new float[2];
		res[0] = (float)Math.cos(phase)*abs;
		res[1] = (float)Math.sin(phase)*abs;
		return res;
	}
}
