

package gps.acquisition;

import cgra.pe.PETrigonometry;

public class AcquisitionPlayground {
	
	/*
	 * Necessary infrmation forthe compuation - rest is left out due to unnecessary overhead
	 */
		private final static int nrFrequencies = 11; //Math.floorDiv((maxRate-minRate), stepRate)+1;
		private final static int stepRate 	= 1000; //1 kHz
		private final static int minRate 		= -5000; //-5kHz
	
	//private final static float[] allFrequencies = {-5000, -4000, -3000, -2000, -1000, 0, 1000, 2000, 3000, 4000, 5000};
	
		/*
		 * Pre-computed multiplications
		 */
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
	
	
	/*
	 * Pre-Computation of frequency shifts
	 */
		private final static float[][] precomputedFrequencyShiftsReal = null;
		private final static float[][] precomputedFrequencyShiftsImag = null;
		private final static int nrPrecomputedFrequencyShifts = -1;
	
		
	/*
	 * Pre-Computed Eulerians for DFTs
	 */
		private final static float[][] precomputedDFTEuleriansReal = null;
		private final static float[][] precomputedDFTEuleriansImag = null;
		
	/*
	 * Find valid acquisition
	 */
		private final static float gamma 		= (float)0.015;
	
	/*
	 * samples and codes
	 */
		private final float[] samplesReal;
		private final float[] samplesImag;
	
		private final float[] codesReal;
		private final float[] codesImag;
	
		private final float[] codesDFTReal;
		private final float[] codesDFTImag;
	
	
	/*
	 *  numbers for samples, FFT
	 */
		private final int nrOfSamples;
		private int nExponent2		= 1;
		private final int N;
		
	/*
	 * FFT Binary reversed indices
	 */
		private final int[] indexPointersBitReversed;
	
	//utlitiy
		private final float oneDivN;
	
	/*
	 * total signal power
	 */
		private float pInTimesNrOfSamples = 0;
	
	/*
	 * Counters for samples
	 */
		private int enteredSamples = 0;
		private int enteredCodes;
	
	/*
	 *  results of the aqcuisition
	 */
		private int dopplerShift = 0;
		private int codeShift    = 0;
	
	
	/**
	 * Constructor with the desired number of codes.
	 * @param nrOfSamples
	 */
	public AcquisitionPlayground(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
		this.oneDivN = 1/nrOfSamples;
		
		enteredCodes = 0;
		
		int nextPowerOfTwo = 2;
		while(nextPowerOfTwo<nrOfSamples) {
			nextPowerOfTwo*=2;
			nExponent2+=1;
		}//find exponent for N
		N = nextPowerOfTwo;
		System.out.println("nrOfSamples:"+nrOfSamples);
		System.out.println("N:"+N);
		System.out.println("exponent:"+nExponent2);
		
		
		/*
		 * Note: automatic zero padding due to Java's zero-initial values in arrays.
		 */
		samplesReal = new float[nrOfSamples];
		samplesImag = new float[nrOfSamples];
		
		codesReal = new float[nrOfSamples];
		codesImag = new float[nrOfSamples];
		
		codesDFTReal = new float[nrOfSamples];
		codesDFTImag = new float[nrOfSamples];
		
		indexPointersBitReversed = new int[nrOfSamples];
		for (int origIndex = 0; origIndex < nrOfSamples; origIndex++) {
			int rev = origIndex;
			int res = 0;
			while (rev!=0){
				  res<<=1;
				  res|=( rev &1);
				  rev>>=1;
				}
			final int revFinal = res;
			indexPointersBitReversed[origIndex] = revFinal;
		}
		
		if (N>nrPrecomputedFrequencyShifts) {
			//TODO handle rest of precomputation here
		}
		
		
	}//Aqcuisition
	
	
	/**
	 * Enters a sample.
	 * Will produce an ArrayIndexOutOfBoundsException if more samples than specified by the constructor are entered.
	 * This exception is left out to avoid overhead
	 * @param real the direct part
	 * @param imag the quadrature part
	 */
	public void enterSample(float real, float imag){
		pInTimesNrOfSamples += (real*real)+(imag*imag);	
		samplesReal[enteredSamples] = real;
		samplesImag[enteredSamples]   = imag;
		enteredSamples++;
	} //enterSample
	
	/**
	 * Enters a code.
	 * Does the transponation of the code vector integrated.
	 * Will produce an ArrayIndexOutOfBoundsException if more codes than specified by the constructor are entered.
	 * This exception is left out to avoid overhead
	 * @param real the direct part
	 * @param imag the quadrature part
	 */
	public void enterCode(float real, float imag){
		codesReal[enteredCodes] = real;
		codesImag[enteredCodes] = -imag;
		enteredCodes++;
	} //enterCode
	
	public boolean startAcquisition(){
		
		int maxFdIndex = 0;
		int maxTauIndex = 0;
		float maxSignalPowerTimesN = 0;
		
		/*
		 * compute DFT of code vector
		 */
		computeDFTCodesOriginal();
		
		/*
		 * Loop: Do the following for each frequency 
		 */
		for (int fd = 0; fd < nrFrequencies; fd++) {
			
		} //fd
		
		//prepare 
		codeShift = maxTauIndex+1;
		dopplerShift = maxFdIndex*stepRate + minRate; //TODO: check if faster than array access to allFrequencies
		return maxSignalPowerTimesN > pInTimesNrOfSamples*gamma ? true : false;
	}
	
	/**
	 * Computes the DFT of the codes for later use. The computation results are stored via side effects.
	 * I.e. codesDFTReal and codesDFTImag are changed.
	 */
	public void computeCodesDFTIteratively() {
		//1st: apply reordering to array
		for (int index = 0; index < N; index++) {
			int shiftedIndex = indexPointersBitReversed[index];
			codesDFTReal[index] = codesReal[shiftedIndex];
			codesDFTImag[index] = codesImag[shiftedIndex];
		}
		//2nd: update according to butterfly
		for (int ii = 2; ii <=N; ii*=2) { // can be left out, nothing to do in that step
			
			//for kk=0, skip long computation of phase
			
			for (int jj = 0; jj < N/ii; jj++) {
				float tmpReal = codesDFTReal[jj*ii+ii/2];
				float tmpImag = codesDFTImag[jj*ii+ii/2];
				
				float realAcc = tmpReal;
				float imagAcc = tmpImag;
				//final +- assignment
				codesDFTReal[jj*ii + ii/2] = codesDFTReal[jj*ii]-realAcc;
				codesDFTImag[jj*ii + ii/2] = codesDFTImag[jj*ii]-imagAcc;
				
				codesDFTReal[jj*ii] +=realAcc;
				codesDFTImag[jj*ii] +=imagAcc; 
			}
			
			//then, the original loop
			for (int kk = 1; kk < ii/2; kk++) {
				
				//compute phase
				final float phase = -2*kk*(1/ii)*(float)Math.PI; 
				final float phaseReal = PETrigonometry.cos(phase);
				final float phaseImag = PETrigonometry.sin(phase);
				
				for (int jj = 0; jj < N/ii; jj++) {
					float tmpReal = codesDFTReal[jj*ii+kk+ii/2];
					float tmpImag = codesDFTImag[jj*ii+kk+ii/2];
					
					float realAcc = tmpReal*phaseReal - tmpImag*phaseImag;
					float imagAcc = tmpReal*phaseImag + tmpImag*phaseReal;
					//final +- assignment
					codesDFTReal[jj*ii + kk + ii/2] = codesDFTReal[jj*ii+kk]-realAcc;
					codesDFTImag[jj*ii + kk+ ii/2]  = codesDFTImag[jj*ii+kk]-imagAcc;
					
					codesDFTReal[jj*ii + kk] +=realAcc;
					codesDFTImag[jj*ii + kk] +=imagAcc; 
				}
			}
		}
		
		
		// greedy implementation
		/*
		for (int kk = 0; kk < N; kk++) {
			//dft[kk] = new float[] {0, 0};
			for (int nn = 0; nn < N; nn++) {
				final float a = codesReal[nn];
				final float b = codesImag[nn];
				final float c = precomputedDFTEuleriansReal[nn][kk];
				final float d = precomputedDFTEuleriansImag[nn][kk];
				final float accReal = a*c - d*b;
				final float accImag = a*d + c*b;
				codesDFTReal[kk]+= accReal;
				codesDFTImag[kk]+= accImag;
			}
		}
		*/
	}
	
	
	/**
	 * Computes the DFT as mixed FFT on the codes vector, of original length
	 */
	public void computeDFTCodesOriginal() {
		
		
	}//
	
	
	
	/*
	 * Getters 
	 */
	public int getDopplerverschiebung(){
		return dopplerShift;
	} //getDopplerverschiebung
	
	public int getCodeVerschiebung(){
		return codeShift;
	}// getCodeVerschiebung
	

	

}
