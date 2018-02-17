package gps.acquisition;

import cgra.pe.PETrigonometry;

public class AcquisitionNoOpt {
	
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
		// first index: fd; second index: nn
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
	public AcquisitionNoOpt(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
		this.oneDivN = 1/nrOfSamples;
		
		enteredCodes = 0;
		
		int nextPowerOfTwo = 2;
		while(nextPowerOfTwo<nrOfSamples) {
			nextPowerOfTwo*=2;
			nExponent2+=1;
		}//find exponent for N
		
		/*
		 * In order to work with classical fast fourier transform (radix-2 cooley tukey, iterative),
		 * N must be a power of 2.
		 * N must also be chosen to be at least as large as 2*nrOfSamples, so that the result of
		 * ifft(fft(samples) x fft(codes*) is not a distorted result, but is exactly equal to the linear convolution
		 * of samples * codes. This linear convolution may be transformed into the circular (desired) convolution by 
		 * time domain aliasing. 
		 * 
		 * The goal of this approach is to take advantage of the full radix-2 fft for the fast convolution and perform that step in N*log2(n)
		 * time; ultimately, the choice of N = 2^(ceil(log2(nrOfSamples))+1) instead of a standard convolution/dft with (not necessarily power
		 * -of-2) nrOfSamples with O(nrOfSamples^2)
		 */
		N = nextPowerOfTwo*2;
		nExponent2++;
		System.out.println("nrOfSamples:"+nrOfSamples);
		System.out.println("N:"+N);
		System.out.println("exponent:"+nExponent2);
		
		
		/*
		 * Note: automatic zero padding due to Java's zero-initial values in arrays.
		 */
		samplesReal = new float[N];
		samplesImag = new float[N];
		
		codesReal = new float[N];
		codesImag = new float[N];
		
		codesDFTReal = new float[N];
		codesDFTImag = new float[N];
		
		indexPointersBitReversed = new int[N];
		for (int origIndex = 0; origIndex < N; origIndex++) {
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
		
		//1st: apply reordering to array
				for (int index = 0; index < N; index++) {
					int shiftedIndex = indexPointersBitReversed[index];
					codesDFTReal[index] = codesReal[shiftedIndex];
					codesDFTImag[index] = codesImag[shiftedIndex];
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
		int index = indexPointersBitReversed[enteredCodes];
		codesDFTReal[index] = real;
		codesDFTImag[index] = -imag;
		enteredCodes++;
	} //enterCode
	
	public boolean startAcquisition(){
		
		int maxFdIndex = 0;
		int maxTauIndex = 0;
		float maxSignalPowerTimesN = 0;
		
		float[] currentFdVectorReal;
		float[] currentFdVectorImag;
		
		
		
		/*
		 * compute DFT of code vector
		 */
		computeCodesDFTIteratively();
		
		/*
		 * Loop: Do the following for each frequency 
		 */
		for (int fd = 0; fd < nrFrequencies; fd++) {
			
			//inthis iteration, work on the precomputed factors for the samples of frequency[fd]
			currentFdVectorReal = precomputedFrequencyShiftsReal[fd];
			currentFdVectorImag = precomputedFrequencyShiftsImag[fd];
			
			
			
			/*
			 * multiplication with precomputed frequency shifts 
			 */
			for (int nn = 0; nn < nrOfSamples; nn++) {
				float factorReal = currentFdVectorReal[nn];
				// real part
					currentFdVectorReal[nn] *= samplesReal[nn];
					currentFdVectorReal[nn] -= samplesImag[nn]*currentFdVectorImag[nn];
				//imagninary part
					currentFdVectorImag[nn] *= samplesReal[nn];
					currentFdVectorImag[nn] += factorReal*samplesImag[nn];
			}//multiplication
			
			//artificial padding, because its faster than multiplication with zero
			for (int nno = nrOfSamples; nno < N; nno++) {
				//apply zero padding afterwards 
				currentFdVectorReal[nno] = 0;
				currentFdVectorImag[nno] = 0;
			}
	
			
			
			/*
			 *  computation of fft of multiplied samples for current fd frequency
			 */
			
			//butterflying
			//TODO check if manual re-ordering of array (as butterfly) or just modified indexing is better
			for (int nn = 0; nn < N; nn++) {
				currentFdVectorReal[nn] = currentFdVectorReal[indexPointersBitReversed[nn]];
				currentFdVectorImag[nn] = currentFdVectorImag[indexPointersBitReversed[nn]];
			}//butterflying
			
			// update according to butterfly
			for (int ii = 2; ii <=N; ii*=2) { //fft outermost loop: recursion level
				
				//for kk=0, skip long computation of phase
				
				for (int jj = 0; jj < N/ii; jj++) {
					float tmpReal = currentFdVectorReal[jj*ii+ii/2];
					float tmpImag = currentFdVectorImag[jj*ii+ii/2];
					
					float realAcc = tmpReal;
					float imagAcc = tmpImag;
					//final +- assignment
					currentFdVectorReal[jj*ii + ii/2] = currentFdVectorReal[jj*ii]-realAcc;
					currentFdVectorImag[jj*ii + ii/2] = currentFdVectorImag[jj*ii]-imagAcc;
					
					currentFdVectorReal[jj*ii] +=realAcc;
					currentFdVectorImag[jj*ii] +=imagAcc; 
				}//for loop "0"-phase
				
				//then, the original loop
				for (int kk = 1; kk < ii/2; kk++) {
					//compute phase
					final float phase = -2*kk*(1/ii)*(float)Math.PI; 
					final float phaseReal = PETrigonometry.cos(phase);
					final float phaseImag = PETrigonometry.sin(phase);
					
					for (int jj = 0; jj < N/ii; jj++) {
						float tmpReal = currentFdVectorReal[jj*ii + kk + ii/2];
						float tmpImag = currentFdVectorImag[jj*ii + kk + ii/2];
						
						float realAcc = tmpReal*phaseReal - tmpImag*phaseImag;
						float imagAcc = tmpReal*phaseImag + tmpImag*phaseReal;
						//final +- assignment
						currentFdVectorReal[jj*ii + kk + ii/2] = currentFdVectorReal[jj*ii+kk]-realAcc;
						currentFdVectorImag[jj*ii + kk+ ii/2]  = currentFdVectorImag[jj*ii+kk]-imagAcc;
						
						currentFdVectorReal[jj*ii + kk] +=realAcc;
						currentFdVectorImag[jj*ii + kk] +=imagAcc; 
					}// jj loop
				}//kk loop
			}//outer loop fft
			
			/*
			 * multiply in frequency domain
			 */
			for (int nn = 0; nn < N; nn++) {
				float aTmp = currentFdVectorReal[nn];
				
				//real part
				currentFdVectorReal[nn] *= codesDFTReal[nn];
				currentFdVectorReal[nn] -= codesDFTImag[nn]*currentFdVectorImag[nn];
				
				//imaginary part
				currentFdVectorImag[nn] *= codesDFTImag[nn];
				currentFdVectorImag[nn] += codesDFTReal[nn]*aTmp;
			}// element-wise multiplication in frequency domain
			
			/*
			 *  inverse fft
			 */
			
			/*
			 *  apply time domain aliasing to rearrange tau indices correctly (downsampling in frequency domain, but more accurate in time domain in our case)
			 */
			
			float currentRFDTimesNMaximum = 0;
			int   currentMaxTauIndex = 0;
			for (int ii =0; ii <nrOfSamples-1; ii++) {
				currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
				currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
				float tmp = currentFdVectorReal[ii]*currentFdVectorReal[ii]
							+ currentFdVectorImag[ii]*currentFdVectorImag[ii];
				if (tmp>currentRFDTimesNMaximum) {
					currentRFDTimesNMaximum = tmp;
					currentMaxTauIndex = ii;
				}// if abs^2 > currentMaximum
			}// aliasing and local maximum detection //TODO: check if maximum extraction is faster in separate loop
			
			//gloabl maximum detection
			if (currentRFDTimesNMaximum > maxSignalPowerTimesN) {
				maxSignalPowerTimesN = currentRFDTimesNMaximum;
				maxFdIndex = fd;
				maxTauIndex = currentMaxTauIndex;
			}


		} //fd

		
		//prepare 
		codeShift = maxTauIndex+1;
		dopplerShift = maxFdIndex*stepRate + minRate; //TODO: check if faster than array access to allFrequencies
		return maxSignalPowerTimesN > pInTimesNrOfSamples*gamma ? true : false;
	}
	
	/**
	 * Computes the DFT of the codes for later use. The computation results are stored via side effects.
	 * I.e. codesDFTReal and codesDFTImag are changed.
	 * Works in-place and requires that the code samples have yet been sorted into the butterfli'ed arrays codesDFTReal/imag
	 */
	public void computeCodesDFTIteratively() {
		
		// update according to butterfly
		for (int ii = 2; ii <=N; ii*=2) { 
			
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
