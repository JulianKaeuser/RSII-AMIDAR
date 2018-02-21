package gps.acquisition;

import cgra.pe.PETrigonometry;

public class AcquisitionWIP {
	
	//atm just for simulatedAnnealing reasons (composition optimization). 
	public static void main (String[] args) {
		AcquisitionWIP acq = new AcquisitionWIP(400);
		acq.startAcquisition();
	}
	
	/*
	 * Necessary infrmation forthe compuation - rest is left out due to unnecessary overhead
	 */
	
		private final static int nrFrequencies = 11; //Math.floorDiv((maxRate-minRate), stepRate)+1;
		private final static int stepRate 	= 1000; //1 kHz
		private final static int minRate 		= -5000; //-5kHz
		
		private final static float fixPrecomputeEFactorPhase = (float)0.07853981633; // = -2*pi*minRate/fs
	
		// pi as constant
		private final static float PI = (float) 3.14159265359;
		private final static float minusPI = (float) -3.14159265359;
		
		/*
		 * Find valid acquisition
		 */
		private final static float gamma 		= (float)0.015;
		
	/*
	 * samples and codes storage
	 */
			private float[] samplesReal;
			private float[] samplesImag;
		
			//holds unprocessed codes after initialization, processed after acquisition
			private float[] codesDFTReal; 
			private float[] codesDFTImag;
			
		
		/*
		 *  numbers for samples, FFT
		 */
			private final int nrOfSamples;
			private int nExponent2	= 1;
			private final int nExponent2Final;
			private final int N;
			
			
			private final boolean radix4;
		
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
	public AcquisitionWIP(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
	        
		int nextPowerOfTwo = 2;
		while(nextPowerOfTwo<nrOfSamples) {
		    nextPowerOfTwo = nextPowerOfTwo << 1;
		    nExponent2+=1;
		}//find exponent for N
		
		/*
		 * In order to work with classical fast fourier transform (radix-2/4 cooley tukey, iterative),
		 * N must be a power of 2.
		 * N must also be chosen to be at least as large as 2*nrOfSamples, so that the result of
		 * ifft(fft(samples) x fft(codes*) is not a distorted result, but is exactly equal to the linear convolution
		 * of samples * codes. This linear convolution may be transformed into the circular (desired) convolution by 
		 * time domain aliasing. 
		 * 
		 * The goal of this approach is to take advantage of the full cooley-tukey fft for the fast convolution and perform that step in N*log2(n)
		 * time; ultimately, the choice of N = 2^(ceil(log2(nrOfSamples))+1) instead of a standard convolution/dft with (not necessarily power
		 * -of-2) nrOfSamples with O(nrOfSamples^2)
		 */
		N = nextPowerOfTwo << 1;
		nExponent2++;
		nExponent2Final = nExponent2;
		radix4 = ((nExponent2 & 0x1) != 1);

		
		    
		    //TODO comment out for final hand in
			 //System.out.println(this.getClass().getName());
		     System.out.println("nrOfSamples:"+nrOfSamples);
		     System.out.println("N:"+N);
		     System.out.println("exponent:"+nExponent2);
		     System.out.println("radix4 "+ (radix4 ? "possible" : "impossible"));
		    
		    //enteredCodes = N-1;
		    enteredCodes = nrOfSamples-1;
		    //enteredCodes = 0;
		    
		    /*
		     * Note: automatic zero padding due to Java's zero-initial values in arrays.
		     */
		    samplesReal = new float[N];
		    samplesImag = new float[N];

		    codesDFTReal =  new float[N];
		    codesDFTImag =  new float[N];
		    

	}//constructor
	
	
	/**
	 * Enters a sample.
	 * Will produce an ArrayIndexOutOfBoundsException if more samples than specified by the constructor are entered.
	 * This exception is left out to avoid overhead
	 * @param real the direct partss
	 * @param imag the quadrature part
	 */
	public void enterSample(float real, float imag){
		pInTimesNrOfSamples += (real*real)+(imag*imag);	
		samplesReal[enteredSamples]   = real;
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
		codesDFTReal[enteredCodes] = real;
		codesDFTImag[enteredCodes] = imag;
		enteredCodes--;
		//enteredCodes--;
	} //enterCode
	
	public boolean startAcquisition(){
		int maxFdIndex = 0;
	    int maxTauIndex = 0;
	    float maxSignalPowerTimesN = 0;
	    int[] maxTauIndices = new int[nrFrequencies];
	    float[] maxSignalPowersTimesN = new float[nrFrequencies];
	    boolean acqResult = false;
	    
	    float[][] fdReals = new float[nrFrequencies][];
	    float[][] fdImags = new float[nrFrequencies][];
	    
	    
	    
	    /*
	     * nrOfSamples is such that a fast convolution, as applied in that case, may be replaced by a radix-4 version, which yields significant speedup
	     * over standard cyclic convolution and over radix 2 (25% less comlex multiplications)
	     * 
	     * flow of computations:
	     * 1. compute radix-4-decimation-in-frequency fft for codes. Codes are flipped by input method. T
	     * 2. for each frequency do:
	     *  2.1 multiply sample vector into "currentFdVector" with frequency shift
	     *  2.2 compute radix-4-decimation-in-frequency fft for "currentFdVector"
	     *  2.3 element-wise multiply both ffts into "currentFdVector"
	     *  2.4 apply inverse radix-4-decimation-in-time fft to "currentFdVector"
	     *  2.5 apply aliasing by nrOfSamples and division by N to "currentFdVector"
	     *  2.6 find maximum of "currentFdVector, store value, dopplerShift, codeShift into array
	     * 3. from comparison array, find global maximum; store results.
	     */
	    if (radix4) {
	    	
		    
		    
		    /*
		     * compute DFT of code vector
		     */
		    //printArraysOctaveFormat("codes unprocessed", true, codesDFTReal, codesDFTImag);
		   computeCodesDFTIterativelyRadix4DIF(); //correct so far
		   //printArraysOctaveFormat("after 1 time dif", true, codesDFTReal, codesDFTImag);
		 //  conjCodes(1);
		   
		   //printArraysOctaveFormat("before dit",true, codesDFTReal, codesDFTImag);
		  // computeCodesDFTIterativelyRadix4DIT();
		  // conjCodes(N);
		   
		  // printArraysOctaveFormat("after dif->conj->dit->conj->/=N(",false, codesDFTReal, codesDFTImag);

		//   conjCodes(N);
		//  compareCodes();
		//  System.exit(0);
		   //TODO delete this stuff all when done / before hand in
		   
		   
		   
		     
		    /*
		     * Loop: Do the following for each frequency 
		     */
		   for (int fd = 0; fd < nrFrequencies; fd++) {
			   //System.out.print("omit fd loop");
			   float[] currentFdVectorReal = new float[N];
			   float[] currentFdVectorImag = new float[N];
			   fdReals[fd] = currentFdVectorReal;
			   fdImags[fd] = currentFdVectorImag;
					  
			   
			   
			    //final int frequency = minRate+fd*stepRate;
		    	final float fdEFactorPhase = (float)-0.01570796326*fd; // = (-2*PI*stepRate*1/fs) * fd
		    	final float fdPrecomputeFactorPhase = fdEFactorPhase + fixPrecomputeEFactorPhase;
		    	
		    	/*
		    	 * multiplication with precomputed frequency shifts 
		    	 */
		    	
		    	for (int nn = 0; nn < nrOfSamples; nn++) {
		    		
		    		//multiply phase with nn here; is faster than incrementing by its initial value every loop -> independent loop iterations
		    		float nnPhasePrecompute = fdPrecomputeFactorPhase*nn;
		    		
		    		//evaluate imaginary e function
		    		float factorReal = PETrigonometry.cos(nnPhasePrecompute);
		    		currentFdVectorReal[nn] = factorReal;
		    		
		    		float factorImag = PETrigonometry.sin(nnPhasePrecompute);
		    		currentFdVectorImag[nn] = factorImag;
		    		
		    		//actual multiplication
		    		// real part
		    			currentFdVectorReal[nn] *= samplesReal[nn];
		    			currentFdVectorReal[nn] -= samplesImag[nn]*factorImag;
		    		//imagninary part
		    			currentFdVectorImag[nn] *= samplesReal[nn];
		    			currentFdVectorImag[nn] += factorReal*samplesImag[nn];
		    	}//multiplication
		    	
		    	printArraysOctaveFormat(fd + "= fd ", true, currentFdVectorReal, currentFdVectorImag);
			   
		    	/*
		    	 * Radix 4 FFT of Samples for current Frequency, DIF
		    	 */
		    	         /*
				          * As DIF produces bit reversed index -ordered output from normal ordered input, and vice versa,
			              * the output is not(!) in normal order. This is, in our case, not necessary, if the inverse FFT is a DIT
			              * approach, which would take bit-reverse index ordered input and produces natural ordered output.
				          */
				final int stages = nExponent2Final>>1; //log(2)/2
				final float phaseFactorBaseDIF = 2*minusPI/N;
				int M = N>>2; // M = N/4
				
				
				//outermost loop: recursion level
				for (int pp = 0; pp<stages; pp++) {			
					final int fourPowerLevel = 1 << (pp<<1);
			
					final float phaseFactorP1Base = 2*fourPowerLevel*phaseFactorBaseDIF;
					final float phaseFactorP2Base = fourPowerLevel*phaseFactorBaseDIF;
					final float phaseFactorP3Base = 3*fourPowerLevel*phaseFactorBaseDIF;
					
					// middle loop: index of block
					final int inc = N>>(pp<<1);
					for (int index = 0; index < N; index+=inc) {
						
						//innermost loop: perform on block
						for (int nn = 0; nn < M; nn++) {
							
							final int indexM = index+M; // testing if faster if in outer loop
							final int index2M = index+2*M; //see line above
							final int index3M = index+3*M; //see line above
							
							float resultReal0, resultImag0;
							float resultReal1, resultImag1;
							float resultReal2, resultImag2;
							float resultReal3, resultImag3;
							
							//fetch required samples from memory
							final float b0Real = currentFdVectorReal[nn+index];
							final float b0Imag = currentFdVectorImag[nn+index];
							final float b1Real = currentFdVectorReal[nn+indexM];
							final float b1Imag = currentFdVectorImag[nn+indexM];
							final float b2Real = currentFdVectorReal[nn+index2M];
							final float b2Imag = currentFdVectorImag[nn+index2M];
							final float b3Real = currentFdVectorReal[nn+index3M];
							final float b3Imag = currentFdVectorImag[nn+index3M];
							
							//put together result values
								
								//result0
								resultReal0 = b0Real+b1Real+b2Real+b3Real;
								resultImag0 = b0Imag+b1Imag+b2Imag+b3Imag;
								
								//result1
								resultReal1 = b0Real-b1Real+b2Real-b3Real;
								resultImag1 = b0Imag-b1Imag+b2Imag-b3Imag;
									//multiplication with phase factor
								final float phaseFactor1 = phaseFactorP1Base*nn;
								final float phaseFactor1Real = PETrigonometry.cos(phaseFactor1);
								final float phaseFactor1Imag = PETrigonometry.sin(phaseFactor1);
									
								final float tmp1 = resultReal1;
									resultReal1 *=phaseFactor1Real;
									resultReal1 -=phaseFactor1Imag*resultImag1;
									
									resultImag1 *=phaseFactor1Real;
									resultImag1 +=phaseFactor1Imag*tmp1;
									
								
								//result2
								resultReal2 = b0Real+b1Imag-b2Real-b3Imag;
								resultImag2 = b0Imag-b1Real-b2Imag+b3Real;
									//multiplication with phase factor
								final float phaseFactor2 = phaseFactorP2Base*nn;
								final float phaseFactor2Real = PETrigonometry.cos(phaseFactor2);
								final float phaseFactor2Imag = PETrigonometry.sin(phaseFactor2);
									
								final float tmp2 = resultReal2;
									resultReal2 *=phaseFactor2Real;
									resultReal2 -=phaseFactor2Imag*resultImag2;
									
									resultImag2 *=phaseFactor2Real;
									resultImag2 +=phaseFactor2Imag*tmp2;
								
								//result3
								resultReal3 = b0Real-b1Imag-b2Real+b3Imag;
								resultImag3 = b0Imag+b1Real-b2Imag-b3Real;
									//mulitplication with phase factor
								final float phaseFactor3 = phaseFactorP3Base*nn;
								final float phaseFactor3Real = PETrigonometry.cos(phaseFactor3);
								final float phaseFactor3Imag = PETrigonometry.sin(phaseFactor3);
									
								final float tmp3 = resultReal3;
									resultReal3 *=phaseFactor3Real;
									resultReal3 -=phaseFactor3Imag*resultImag3;
									
									resultImag3 *=phaseFactor3Real;
									resultImag3 +=phaseFactor3Imag*tmp3;
							
							
							//assignments to memory
							currentFdVectorReal[nn+index] = resultReal0;
							currentFdVectorImag[nn+index] = resultImag0;
							
							currentFdVectorReal[nn+indexM] = resultReal1;
							currentFdVectorImag[nn+indexM] = resultImag1;
							
							currentFdVectorReal[nn+index2M] = resultReal2;
							currentFdVectorImag[nn+index2M] = resultImag2;
							
							currentFdVectorReal[nn+index3M] = resultReal3;
							currentFdVectorImag[nn+index3M] = resultImag3;
						}//nn
					} //index
					M = M>>2;
				}//p = recursionLevel
				
				//so far, error of 10^-3... keep testing
				
		    	printArraysOctaveFormat("fd"+fd+" after DIF", true, currentFdVectorReal, currentFdVectorImag);
		    	/*
		    	 * multiply in frequency domain
		    	 */
		    	
		    	for (int nn = 0; nn < N; nn++) {
		    		float aTmp = currentFdVectorReal[nn];
		    		
		    		//real part
		    		currentFdVectorReal[nn] *= codesDFTReal[nn];
		    		currentFdVectorReal[nn] -= codesDFTImag[nn]*currentFdVectorImag[nn]; 
		    		currentFdVectorReal[nn] /=N;
		    		
		    																			 
		    		
		    		//imaginary part
		    		currentFdVectorImag[nn] *= codesDFTReal[nn];  
		    		currentFdVectorImag[nn] += codesDFTImag[nn]*aTmp;
		    		currentFdVectorImag[nn] *= (float) -1; //conjugate complex -> next forward fft will be reversed
		    		currentFdVectorImag[nn] /= N;
		    	}// element-wise multiplication in frequency domain
		    	
		    	printArraysOctaveFormat("fd"+fd+" after freq domain mult", true, currentFdVectorReal, currentFdVectorImag);
		    	
		    	
		    	/*
		    	 * Inverse Radix 4 FFT for current frequency, DIT
		    	 */
		    	/*
				 * The DIT-Approach usually takes bit-reversed index ordered inputs and produces normal-ordered
				 * outputs. Therefore, make sure that the input is bit-reversed index ordered (e.g. by explicit reordering
				 * or creating it from a DIF-Forward FFT (our case)).
				 */
				//final int stages = nExponent2>>1; // log2(N)/2; yet defined
				
		    	
		    	
				final float phaseFactorBaseDIT = minusPI/N; //no minus because of inversion 
				
				final float phaseFactorBase1 = phaseFactorBaseDIT*(float)2.0;
				final float phaseFactorBase2 = phaseFactorBaseDIT*(float)4.0;
				final float phaseFactorBase3 = phaseFactorBaseDIT*(float)6.0;
				
				for (int pp = 0; pp<stages; pp++) {                                    
					final int fourPowerLevel = 4<<(pp<<1);
					
					
					final int offset = 1<<(pp<<1); 
					final int offset2 = offset<<1; //offset*2;
					final int offset3 = 3*offset; //offset *3;
					
					final int phaseFactorInc = N/fourPowerLevel; //TODO exchange with shifts
					for (int index = 0; index < N; index+=fourPowerLevel) {
						int kk = 0;
						for (int nn = index; nn < index+offset; nn++) {
							
							
							
							//phase factors
							final float phaseFactorKK = phaseFactorBase1*kk; 
							final float phaseFactor2KK = phaseFactorBase2*kk; 
							final float phaseFactor3KK = phaseFactorBase3*kk; 
							
							final float phaseFactor1Real = PETrigonometry.cos(phaseFactorKK);
						    final float phaseFactor1Imag = PETrigonometry.sin(phaseFactorKK);
						    
						    final float phaseFactor2Real = PETrigonometry.cos(phaseFactor2KK); 
						    final float phaseFactor2Imag = PETrigonometry.sin(phaseFactor2KK);
						    
						    final float phaseFactor3Real = PETrigonometry.cos(phaseFactor3KK); 
						    final float phaseFactor3Imag = PETrigonometry.sin(phaseFactor3KK);
						    		
							//fetch required samples from memory
							float b0Real = currentFdVectorReal[nn];			
							float b0Imag = currentFdVectorImag[nn];
							float b1Real = currentFdVectorReal[nn+offset];		
							float b1Imag = currentFdVectorImag[nn+offset];
							float b2Real = currentFdVectorReal[nn+offset2];	
							float b2Imag = currentFdVectorImag[nn+offset2];
							float b3Real = currentFdVectorReal[nn+offset3];	
							float b3Imag = currentFdVectorImag[nn+offset3];
							
							// multiplication with phase factors
							
							float a0Real, a0Imag;
							float a1Real, a1Imag;
							float a2Real, a2Imag;
							float a3Real, a3Imag;
							
							//intermediate values (multiplication with phase factors
							a0Real = b0Real;
							a0Imag = b0Imag;
							
							a1Real = b2Real*phaseFactor1Real - b2Imag*phaseFactor1Imag;
							a1Imag = b2Real*phaseFactor1Imag + b2Imag*phaseFactor1Real;
							
							a2Real = b1Real*phaseFactor2Real - b1Imag*phaseFactor2Imag;
							a2Imag = b1Real*phaseFactor2Imag + b1Imag*phaseFactor2Real;
							
							a3Real = b3Real*phaseFactor3Real - b3Imag*phaseFactor3Imag;
							a3Imag = b3Real*phaseFactor3Imag + b3Imag*phaseFactor3Real;
							
							//assign to memory
							currentFdVectorReal[nn] 		 = a0Real + a2Real + a1Real + a3Real;
							currentFdVectorImag[nn] 		 = a0Imag + a2Imag + a1Imag + a3Imag;
														
							currentFdVectorReal[nn+offset]  = a0Real - a2Real + a1Imag - a3Imag;
							currentFdVectorImag[nn+offset]  = a0Imag - a2Imag - a1Real + a3Real;
														
							currentFdVectorReal[nn+offset2] = a0Real + a2Real - a1Real - a3Real; 
							currentFdVectorImag[nn+offset2] = a0Imag + a2Imag - a1Imag - a3Imag;
														
							currentFdVectorReal[nn+offset3] = a0Real - a2Real - a1Imag + a3Imag;
							currentFdVectorImag[nn+offset3] = a0Imag - a2Imag + a1Real - a3Real;
					
							kk +=phaseFactorInc; //increment phase factor rotation
						}//nn loop
					}//index
					//if (fourPowerLevel==N) {
					//	break;
					//}
				}//pp = recursion level	
				
		    	printArraysOctaveFormat("fd"+fd+" after ifft DIT", true, currentFdVectorReal, currentFdVectorImag);


			   /*
		    	 *  apply time domain aliasing to rearrange tau indices correctly 
		    	 *  (downsampling in frequency domain, but more accurate in time domain in our case)
		    	 *       apply scaling by 1/N after inverse FFT here, also
		    	 */
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	int maxTauIndexFd = 0;
		    	int ii;
		    	for (ii = 0; ii <nrOfSamples; ii++) { //TODO check how to implement maximum check for last step 
		    		//aliasing: add up what should have been wrapped
		    		currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
		    		currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
		    		currentFdVectorReal[ii] /= N;
		    		currentFdVectorImag[ii] /= N;
		    		
		    		float tmp  = currentFdVectorReal[ii]*currentFdVectorReal[ii];
		    		float tmp2 = currentFdVectorImag[ii]*currentFdVectorImag[ii];
		    		//tmp  *= nFactorRadix4Squared;
		    		//tmp2 *= nFactorRadix4Squared;
		    		
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxTauIndexFd = ii;
		    			currentRFDTimesNMaximum = tmp;
		    			
		    		}// if abs^2 > currentMaximum
		    		//ii++;
		    	}// aliasing and local maximum detection
		    	
	    		maxTauIndices[fd] = maxTauIndexFd;
	    		maxSignalPowersTimesN[fd] = currentRFDTimesNMaximum;
	    		
	    		printArraysOctaveFormat("fd"+fd+" after aliasing", true, currentFdVectorReal, currentFdVectorImag);
		   }//fd
	    } //radix 4
	    
	    for (int fd = 0; fd < nrFrequencies; fd++) {
	    	float[] absVector = new float[N];
    		float[] nullVector = new float[N];
    		
    		for (int nn = 0; nn < N; nn++) {
    			absVector[nn] = fdReals[fd][nn]*fdReals[fd][nn] + fdImags[fd][nn]*fdImags[fd][nn];
    			absVector[nn] /= (N*N);
    		}
    		printArraysOctaveFormat("a"+fd, true, absVector);
    		
	    	
	    }
		
		
		/*
		 * Equal for both variants:
	     * Find absoulte maximum from all frequencies
	     */
	    for (int fd = 0; fd < nrFrequencies; fd++) {
	    	if (maxSignalPowersTimesN[fd] > maxSignalPowerTimesN) {
	    		maxSignalPowerTimesN = maxSignalPowersTimesN[fd];
	    		maxFdIndex = fd;
	    		maxTauIndex = maxTauIndices[fd];
	    	}
	    }// absolute maximum, factored out of other loop to enable more parallelization
	    
	    // assigments
	    codeShift = maxTauIndex+1;
	    System.out.println("maxFDIndex =" +maxFdIndex);
	    dopplerShift = maxFdIndex*stepRate + minRate; 
	    acqResult =  maxSignalPowerTimesN > pInTimesNrOfSamples*gamma*N*N; //division by N had been omitted in inverse FFT; this way, only needed once
		return acqResult;
	}//startAcquisition
	
	/**
	 * Computes the Radix-4 Decimation-in-Frequency FFT in an iterative way for the codes vector.
	 * In-place method - the result is stored in the code vector itself.
	 * As a prerequisite, the length of the code vector must be a even power of two (= power of four).
	 
	 */
	public void computeCodesDFTIterativelyRadix4DIF() {
		//printArraysOctaveFormat("raw codes", codesDFTReal, codesDFTImag);
		
		/*
		 * As DIF produces bit reversed index -ordered output from normal ordered input, and vice versa,
	     * the output is not(!) in normal order. This is, in our case, not necessary, if the inverse FFT is a DIT
	     * approach, which would take bit-reverse index ordered input and produces natural ordered output.
		 */
		final int stages = nExponent2>>1; //log(2)/2
		final float phaseFactorBase = -2*PI/N;
		int M = N>>2; // M = N/4
		
		//outermost loop: recursion level
		for (int pp = 0; pp<stages; pp++) {			
			final int fourPowerLevel = 1 << (pp<<1);
	
			final float phaseFactorP1Base = 2*fourPowerLevel*phaseFactorBase;
			final float phaseFactorP2Base = fourPowerLevel*phaseFactorBase;
			final float phaseFactorP3Base = 3*fourPowerLevel*phaseFactorBase;
			
			// middle loop: index of block
			final int inc = N>>(pp<<1);
			for (int index = 0; index < N; index+=inc) {
				
				//innermost loop: perform on block
				for (int nn = 0; nn < M; nn++) {
					
					final int indexM = index+M; // testing if faster if in outer loop
					final int index2M = index+2*M; //see line above
					final int index3M = index+3*M; //see line above
					
					float resultReal0, resultImag0;
					float resultReal1, resultImag1;
					float resultReal2, resultImag2;
					float resultReal3, resultImag3;
					
					//fetch required samples from memory
					float b0Real = codesDFTReal[nn+index];
					float b0Imag = codesDFTImag[nn+index];
					float b1Real = codesDFTReal[nn+indexM];
					float b1Imag = codesDFTImag[nn+indexM];
					float b2Real = codesDFTReal[nn+index2M];
					float b2Imag = codesDFTImag[nn+index2M];
					float b3Real = codesDFTReal[nn+index3M];
					float b3Imag = codesDFTImag[nn+index3M];
					
					//put together result values
						
						//result0
						resultReal0 = b0Real+b1Real+b2Real+b3Real;
						resultImag0 = b0Imag+b1Imag+b2Imag+b3Imag;
						
						//result1
						resultReal1 = b0Real-b1Real+b2Real-b3Real;
						resultImag1 = b0Imag-b1Imag+b2Imag-b3Imag;
							//multiplication with phase factor
							float phaseFactor1 = phaseFactorP1Base*nn;
							float phaseFactor1Real = PETrigonometry.cos(phaseFactor1);
							float phaseFactor1Imag = PETrigonometry.sin(phaseFactor1);
							
							float tmp1 = resultReal1;
							resultReal1 *=phaseFactor1Real;
							resultReal1 -=phaseFactor1Imag*resultImag1;
							
							resultImag1 *=phaseFactor1Real;
							resultImag1 +=phaseFactor1Imag*tmp1;
							
						
						//result2
						resultReal2 = b0Real+b1Imag-b2Real-b3Imag;
						resultImag2 = b0Imag-b1Real-b2Imag+b3Real;
							//multiplication with phase factor
							float phaseFactor2 = phaseFactorP2Base*nn;
							float phaseFactor2Real = PETrigonometry.cos(phaseFactor2);
							float phaseFactor2Imag = PETrigonometry.sin(phaseFactor2);
							
							float tmp2 = resultReal2;
							resultReal2 *=phaseFactor2Real;
							resultReal2 -=phaseFactor2Imag*resultImag2;
							
							resultImag2 *=phaseFactor2Real;
							resultImag2 +=phaseFactor2Imag*tmp2;
						
						//result3
						resultReal3 = b0Real-b1Imag-b2Real+b3Imag;
						resultImag3 = b0Imag+b1Real-b2Imag-b3Real;
							//mulitplication with phase factor
							float phaseFactor3 = phaseFactorP3Base*nn;
							float phaseFactor3Real = PETrigonometry.cos(phaseFactor3);
							float phaseFactor3Imag = PETrigonometry.sin(phaseFactor3);
							
							float tmp3 = resultReal3;
							resultReal3 *=phaseFactor3Real;
							resultReal3 -=phaseFactor3Imag*resultImag3;
							
							resultImag3 *=phaseFactor3Real;
							resultImag3 +=phaseFactor3Imag*tmp3;
					
					
					//assignments to memory
					codesDFTReal[nn+index] = resultReal0;
					codesDFTImag[nn+index] = resultImag0;
					
					codesDFTReal[nn+indexM] = resultReal1;
					codesDFTImag[nn+indexM] = resultImag1;
					
					codesDFTReal[nn+index2M] = resultReal2;
					codesDFTImag[nn+index2M] = resultImag2;
					
					codesDFTReal[nn+index3M] = resultReal3;
					codesDFTImag[nn+index3M] = resultImag3;
				}//nn
			} //index
			M = M>>2;
		}//p = recursionLevel
		
		//printArraysOctaveFormat("dif-fft of codes:", codesDFTReal, codesDFTImag);
	}//computeCodesDFTIterativelyRadix4

	public void conjCodes(float divFactor) {
		for (int nn = 0; nn < N; nn++) {
			codesDFTImag[nn] *= -1;
			codesDFTReal[nn] /=divFactor;
			codesDFTImag[nn] /= divFactor;
		}
	}
	
	public static void printArraysOctaveFormat(String msg, boolean octave, float[] real, float[] imag) {
		if (msg.length() !=0) return;
		if(octave) {
			System.out.println(msg);
			String st = " = [";
				for (int nn = 0; nn < real.length; nn++) {
					st+= "("+real[nn] +" + "+imag[nn]+ "i)	 ";
				}
				st += "];";
				System.out.println(st);
		}
		
		else {
			System.out.println(msg);
			String st = "";
			for (int nn = 0; nn < real.length; nn++) {
				st= nn + " ("+real[nn] +" + "+imag[nn]+ "i)	 ";
				System.out.println(st);

			}
			
		}
	}
	
public static void printArraysOctaveFormat(String msg, boolean octave, float[] abs) {
	if (msg.length() !=0) return;
		
		
		if(octave) {
			//System.out.print(msg);
			String st = msg+ "= [";
				for (int nn = 0; nn < abs.length; nn++) {
					st+= abs[nn]+ " ";
				}
				st += "];";
				System.out.println(st);
		}
		
		else {
			System.out.println(msg);
			String st = "";
			for (int nn = 0; nn < abs.length; nn++) {
				st= nn + " " + abs[nn];
				System.out.println(st);

			}
			
		}
	}
	
	/**
	 * For Testing the inverse DIT approach
	 * Computes the inverse FFT (without scaling 1/N) of the codes vector with an iterative Decimation-in-Time radix-4
	 * approach. 
	 */
	public void computeCodesDFTIterativelyRadix4DIT() {
		
		
		/*
		 * The DIT-Approach usually takes bit-reversed index ordered inputs and produces normal-ordered
		 * outputs. Therefore, make sure that the input is bit-reversed index ordered (e.g. by explicit reordering
		 * or creating it from a DIF-Forward FFT (our case)).
		 */
		final int stages = nExponent2>>1; // log2(N)/2;
		
		final float phaseFactorBase = minusPI/N; //no minus because of inversion 
		
		final float phaseFactorBase1 = phaseFactorBase*(float)2.0;
		final float phaseFactorBase2 = phaseFactorBase*(float)4.0;
		final float phaseFactorBase3 = phaseFactorBase*(float)6.0;
		
		for (int pp = 0; pp<stages; pp++) {                                    
			final int fourPowerLevel = 4<<(pp<<1);
			
			
			final int offset = 1<<(pp<<1); 
			final int offset2 = offset<<1; //offset*2;
			final int offset3 = 3*offset; //offset *3;
			
			final int phaseFactorInc = N/fourPowerLevel; //TODO exchange with shifts
			for (int index = 0; index < N; index+=fourPowerLevel) {
				int kk = 0;
				for (int nn = index; nn < index+offset; nn++) {
					
					
					
					//phase factors
					final float phaseFactorKK = phaseFactorBase1*kk; 
					final float phaseFactor2KK = phaseFactorBase2*kk; 
					final float phaseFactor3KK = phaseFactorBase3*kk; 
					
					final float phaseFactor1Real = PETrigonometry.cos(phaseFactorKK);
				    final float phaseFactor1Imag = PETrigonometry.sin(phaseFactorKK);
				    
				    final float phaseFactor2Real = PETrigonometry.cos(phaseFactor2KK); 
				    final float phaseFactor2Imag = PETrigonometry.sin(phaseFactor2KK);
				    
				    final float phaseFactor3Real = PETrigonometry.cos(phaseFactor3KK); 
				    final float phaseFactor3Imag = PETrigonometry.sin(phaseFactor3KK);
				    		
					//fetch required samples from memory
					float b0Real = codesDFTReal[nn];			
					float b0Imag = codesDFTImag[nn];
					float b1Real = codesDFTReal[nn+offset];		
					float b1Imag = codesDFTImag[nn+offset];
					float b2Real = codesDFTReal[nn+offset2];	
					float b2Imag = codesDFTImag[nn+offset2];
					float b3Real = codesDFTReal[nn+offset3];	
					float b3Imag = codesDFTImag[nn+offset3];
					
					// multiplication with phase factors
					
					float a0Real, a0Imag;
					float a1Real, a1Imag;
					float a2Real, a2Imag;
					float a3Real, a3Imag;
					
					//intermediate values (multiplication with phase factors)
					a0Real = b0Real;
					a0Imag = b0Imag;
					
					
					a1Real = b2Real*phaseFactor1Real - b2Imag*phaseFactor1Imag;
					a1Imag = b2Real*phaseFactor1Imag + b2Imag*phaseFactor1Real;
					
					
					a2Real = b1Real*phaseFactor2Real - b1Imag*phaseFactor2Imag;
					a2Imag = b1Real*phaseFactor2Imag + b1Imag*phaseFactor2Real;
					
					
					a3Real = b3Real*phaseFactor3Real - b3Imag*phaseFactor3Imag;
					a3Imag = b3Real*phaseFactor3Imag + b3Imag*phaseFactor3Real;
					
					//assign to memory
					codesDFTReal[nn] 		 = a0Real + a2Real + a1Real + a3Real;
					codesDFTImag[nn] 		 = a0Imag + a2Imag + a1Imag + a3Imag;
					
					codesDFTReal[nn+offset]  = a0Real - a2Real + a1Imag - a3Imag;
					codesDFTImag[nn+offset]  = a0Imag - a2Imag - a1Real + a3Real;
					
					codesDFTReal[nn+offset2] = a0Real + a2Real - a1Real - a3Real; 
					codesDFTImag[nn+offset2] = a0Imag + a2Imag - a1Imag - a3Imag;
					
					codesDFTReal[nn+offset3] = a0Real - a2Real - a1Imag + a3Imag;
					codesDFTImag[nn+offset3] = a0Imag - a2Imag + a1Real - a3Real;
			
				    		
					kk +=phaseFactorInc;
				}//nn loop
			}//index
			//if (fourPowerLevel==N) {
			//	break;
			//}
		}//pp = recursion level	
		
		//printArraysOctaveFormat("codes after DIT", codesDFTReal, codesDFTImag);
	}//computeCodesDFTIterativelyRadix4DIT
	   
	public void compareCodes() {
		float byN = 1/N;
		
		//aliasing ? 
		//for (int nn=0; nn < nrOfSamples-1; nn++) {
		//	codesDFTReal[nn] +=codesDFTReal[nn+nrOfSamples];
		//	codesDFTImag[nn] +=codesDFTImag[nn+nrOfSamples];
		//}
		//for (int nn = 0; nn < N; nn++) {
		//	float realDiff = codesDFTReal[nn] - codesRefReal[nn];
		//	float imagDiff = codesDFTImag[nn] - codesRefImag[nn];
		//	String line = nn+ ":\tref: ("+codesRefReal[nn]+","+codesRefImag[nn]+")  \t computed:("+codesDFTReal[nn]+","+codesDFTImag[nn]+")\t diff: ("+realDiff+","+imagDiff+")";
		//	//String line = codesDFTReal[nn] + " + "+ codesDFTImag[nn]+"i";
		//	System.out.println(line);
		//}
	}
	
	public void computeCodesDFTIteratively() {
		
		// update according to butterfly
		int p = 1;
		for (int ii = 2; ii <=N; ii*=2) { 
			final int halfII = ii<<1;
			final float phaseConst = -2*(1/ii)*PI; 
			//for kk=0, skip long computation of phase
			final int nByII = N>>p;
			p++;
			
			
			//then, the original loop
			final int iiBy2 = ii<<1;
			for (int kk = 0; kk < iiBy2; kk++) {
				
				//compute phase
				final float phase = phaseConst*kk;
				final float phaseReal = PETrigonometry.cos(phase);
				final float phaseImag = PETrigonometry.sin(phase);
				
				for (int jj = 0; jj < nByII; jj++) {
					float tmpReal = codesDFTReal[jj*ii+kk+halfII];
					float tmpImag = codesDFTImag[jj*ii+kk+halfII];
					
					float realAcc1 = tmpReal*phaseReal;
					float imagAcc1 = tmpReal*phaseImag;
					float realAcc2 = tmpImag*phaseImag;
					float imagAcc2 = tmpImag*phaseReal;
					//final +- assignment
					float aReal = codesDFTReal[jj*ii+kk]-realAcc1; 
				    float aImag = codesDFTImag[jj*ii+kk]-imagAcc1; 
					aReal +=realAcc2;
					aImag -= imagAcc2;
					codesDFTReal[jj*ii + kk + halfII] = aReal;
					codesDFTImag[jj*ii + kk+ halfII]  = aImag;
					
					codesDFTReal[jj*ii + kk] +=realAcc1;
					codesDFTImag[jj*ii + kk] +=imagAcc1; 
					codesDFTReal[jj*ii + kk] -=realAcc2;
					codesDFTImag[jj*ii + kk] +=imagAcc2; 
				}
			}
		}
	}	//codesDFT

	
	
	
	
/*
 * Getters for the output data
 */

	/**
	 * Returns the computed doppler shift
	 * @return the doppler shift
	 */
	
	public int getDopplerverschiebung(){
		return dopplerShift;
	}//getDopplerverschiebung

	/**
	 * Returns the found code shift
	 * @return the code shift
	 */
	public int getCodeVerschiebung(){
		return codeShift;
	} // getCodeVerschiebung

}//end of class
