package gps.acquisition;


public class Acquisition {
	
	private final static float samplingRate = 400000; //400kHz
	
	private final static int nrFrequencies = 11; //Math.floorDiv((maxRate-minRate), stepRate)+1;
	
	private final static float[] allFrequencies = {-5000, -4000, -3000, -2000, -1000, 0, 1000, 2000, 3000, 4000, 5000};
	
	private final static float[] allFrequenciesDivSamplingRateTimes2Pi =   {(float)-0.07853982, 
																			(float)-0.06283186, 
																			(float)-0.04712389, 
																			(float)-0.03141593, 
																			(float)-0.015707964, 
																			(float)0.0, 
																			(float)0.015707964, 
																			(float)0.03141593, 
																			(float)0.04712389, 
																			(float)0.06283186, 
																			(float)0.07853982
																			};
	
	private final static float[][] precomputedFrequencyShiftsReal = null;
	private final static float[][] precomputedFrequencyShiftsImag = null;
	private final static int nrPrecomputedFrequencies = -1;
	
	private final static float gamma 		= (float)0.015;
	
	//samples and codes
	private final float[] samplesReal;
	private final float[] samplesImag;
	
	
	// numbers for samples, FFT
	private final int nrOfSamples;
	private int nExponent2		= 1;
	private int N 				= 1;
	
	//utlitiy
	private final float oneDivNrOfSamples;
	private float oneDivN = 1;
	
	//total signal power
	private float pInTimesNrOfSamples = 0;
	
	public Acquisition(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
		oneDivN = 1/nrOfSamples;
		
		N = 1;
		
		while(nrOfSamples>N) {
			N*=2;
			nExponent2 +=1;
		}
		
		samplesReal = new float[N];
		samplesImag = new float[N];
		
		
	}
	
	public void enterSample(float real, float imag){
		pInTimesNrOfSamples += (real*real)+(imag*imag);	
	}
	
	public void enterCode(float real, float imag){
	}
	
	public boolean startAcquisition(){
		
		return false;
	}
	
	public int getDopplerverschiebung(){
		return 0;
	}
	
	public int getCodeVerschiebung(){
		return 0;
	}
	
	

}
