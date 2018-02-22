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
		private final static float minus2PI = (float) 2*minusPI;
		
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
	} //enterCode
	
	public boolean startAcquisition(){
		int maxFdIndex = 0;
	    int maxTauIndex = 0;
	    float maxSignalPowerTimesN = 0;
	    int[] maxTauIndices = new int[nrFrequencies];
	    float[] maxSignalPowersTimesN = new float[nrFrequencies];
	    boolean acqResult = false;
	    

	    
	    
	    
	    /*
	     * nrOfSamples is such that a fast convolution, as applied in that case, may be replaced by a radix-4 version, which yields significant speedup
	     * over standard cyclic convolution and over radix 2 (25% less comlex multiplications)
	     * 
	     * flow of computations:
	     * 1. compute radix-4-decimation-in-frequency fft for codes. Codes are flipped by input method. T
	     * 2. for each frequency do:
	     *  2.1 multiply sample vector into "currentFdVector" with frequency shift
	     *  2.2 compute radix-4-decimation-in-frequency fft for "currentFdVector"
	     *  2.3 element-wise multiply both ffts into "currentFdVector", complex conjugate the result
	     *  2.4 apply inverse radix-4-decimation-in-time fft to "currentFdVector"
	     *  2.5 apply aliasing by nrOfSamples and division by N to "currentFdVector"
	     *  2.6 find maximum of "currentFdVector, store value, dopplerShift, codeShift into array
	     * 3. from comparison array, find global maximum; store results.
	     */
	    if (radix4) {
		    
		    /*
		     * compute DFT of code vector
		     */
		   computeCodesDFTIterativelyRadix4DIF(); //correct so far		    
		     
		    /*
		     * Loop: Do the following for each frequency 
		     */
		   for (int fd = 0; fd < nrFrequencies; fd++) {
			   float[] currentFdVectorReal = new float[N];
			   float[] currentFdVectorImag = new float[N];
					  
		    	final float fdEFactorPhase = (float)-0.01570796326*fd; // = (-2*PI*stepRate*1/fs) * fd
		    	final float fdPrecomputeFactorPhase = fdEFactorPhase + fixPrecomputeEFactorPhase;
		    	
		    	/*
		    	 * multiplication with precomputed frequency shifts 
		    	 */
		    	
		    	for (int nn = 0; nn < nrOfSamples; nn++) {
		    		
		    		//multiply phase with nn here;
		    		//is faster than incrementing by its initial value every loop -> independent loop iterations
		    		// should also be better for precision 
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
		    	
			   
		    	/*
		    	 * Radix 4 FFT of Samples for current Frequency, DIF
		    	 */
		    	         /*
				          * As DIF produces bit reversed index -ordered output from normal ordered input, and vice versa,
			              * the output is not(!) in normal order. This is, in our case, not necessary, if the inverse FFT is a DIT
			              * approach, which would take bit-reverse index ordered input and produces natural ordered output.
				          */
				final int stages = nExponent2Final>>1; //log(2)/2
				final float phaseFactorBaseDIF = minus2PI/N;
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
				
				
		    	/*
		    	 * multiply in frequency domain
		    	 */
		    	
		    	for (int nn = 0; nn < N; nn++) {
		    		float aTmp = currentFdVectorReal[nn];
		    		
		    		//real part
		    		currentFdVectorReal[nn] *= codesDFTReal[nn];
		    		currentFdVectorReal[nn] -= codesDFTImag[nn]*currentFdVectorImag[nn]; 
 		
		    		//imaginary part
		    		currentFdVectorImag[nn] *= codesDFTReal[nn];  
		    		currentFdVectorImag[nn] += codesDFTImag[nn]*aTmp;
		    		currentFdVectorImag[nn] *= -1;  //apply complex conjugation to use forward fft as inverse
		    	}// element-wise multiplication in frequency domain
		    	
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
					
					final int phaseFactorInc = (N>>((pp<<1)+2)); //TODO exchange with shifts
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
				
		    	


			    /*
		    	 *  apply time domain aliasing to rearrange tau indices correctly 
		    	 *  (downsampling in frequency domain, but more accurate in time domain in our case)
		    	 */
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	int maxTauIndexFd = 0;
		    	int ii;
		    	for (ii = 0; ii <nrOfSamples; ii++) { 
		    		//aliasing: add up what should have been wrapped
		    		currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
		    		currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
		    		currentFdVectorReal[ii] /= N;
		    		currentFdVectorImag[ii] /= N;
		    		
		    		float tmp  = currentFdVectorReal[ii]*currentFdVectorReal[ii];
		    		float tmp2 = currentFdVectorImag[ii]*currentFdVectorImag[ii];
		    		
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxTauIndexFd = ii;
		    			currentRFDTimesNMaximum = tmp;
		    			
		    		}// if abs^2 > currentMaximum
		    		//ii++;
		    	}// aliasing and local maximum detection
		    	
	    		maxTauIndices[fd] = maxTauIndexFd;
	    		maxSignalPowersTimesN[fd] = currentRFDTimesNMaximum;
		   }//fd
	    } //radix 4
	    
	    /*
	     * nrOfSamples is such that a fast convolution, as applied in that case, must only be replaced by a radix-2 version,  
	     * flow of computations:
	     * 1. compute radix-2-decimation-in-frequency fft for codes. Codes are flipped by input method. 
	     * 2. for each frequency do:
	     *  2.1 multiply sample vector into "currentFdVector" with frequency shift
	     *  2.2 compute radix-2-decimation-in-frequency fft for "currentFdVector"
	     *  2.3 element-wise multiply both ffts into "currentFdVector", complex conjugate the result
	     *  2.4 apply inverse radix-2-decimation-in-time fft to "currentFdVector"
	     *  2.5 apply aliasing by nrOfSamples and division by N to "currentFdVector"
	     *  2.6 find maximum of "currentFdVector, store value, dopplerShift, codeShift into array
	     * 3. from comparison array, find global maximum; store results.
	     */
	    else { //radix 2
	    	
	    	/*
	    	 * Compute DFT of Codes vector
	    	 */
	    	computeCodesDFTIterativelyRadix2DIF();
	    	
	    		    	
	    	/*
	    	 * For each frequency do: 
	    	 */
	    	for (int fd = 0; fd < nrFrequencies; fd++) {
	    		float[] currentFdVectorReal = new float[N];
	    		float[] currentFdVectorImag = new float[N];
	    		
	    		// multiplication with unity root(fd)
	    		final float fdEFactorPhase = (float)-0.01570796326*fd; // = (-2*PI*stepRate*1/fs) * fd
		    	final float fdPrecomputeFactorPhase = fdEFactorPhase + fixPrecomputeEFactorPhase;
		    	
		    	/*
		    	 * multiplication with precomputed frequency shifts 
		    	 */
		    	
		    	for (int nn = 0; nn < nrOfSamples; nn++) {
		    		
		    		//multiply phase with nn here;
		    		//is faster than incrementing by its initial value every loop -> independent loop iterations
		    		// should also be better for precision 
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
		    	
		    	/*
		    	 * Radix-2 DIF FFT of sample vector
		    	 */
		    	final int stages = nExponent2Final;
				final float phaseFactorExp = minus2PI/N;
				
				for (int pp = 0; pp < stages; pp++) {
					final int twoPowerLevel = 1<<(pp);
					final int halfsize = N>>(pp+1);
					final int size = N>>(pp);
					
					for (int index = 0; index < twoPowerLevel; index++) {
						final int beginKK = index*size;
						final int endKK = beginKK+halfsize;
						int phaseFactorInc = 0;
						for (int nn = beginKK; nn < endKK; nn++) {
							 
							
							final float tmpReal = currentFdVectorReal[nn];
							final float tmpImag = currentFdVectorImag[nn];
							
							final float phaseFactor = phaseFactorInc*phaseFactorExp;
							final float phaseFactorReal = PETrigonometry.cos(phaseFactor);
							final float phaseFactorImag = PETrigonometry.sin(phaseFactor);
							
							final float a0Real = tmpReal - currentFdVectorReal[nn+halfsize];
							final float a0Imag = tmpImag - currentFdVectorImag[nn+halfsize];
							
							float b0Real  = a0Real*phaseFactorReal;
								  b0Real -= a0Imag*phaseFactorImag;
								  
							float b0Imag  = a0Imag*phaseFactorReal;
							      b0Imag += a0Real*phaseFactorImag;
							
							//assignments to memory
							currentFdVectorReal[nn] = tmpReal + currentFdVectorReal[nn+halfsize];
							currentFdVectorImag[nn] = tmpImag + currentFdVectorImag[nn+halfsize];
							
							currentFdVectorReal[nn+halfsize] = b0Real;
							currentFdVectorImag[nn+halfsize] = b0Imag;
							
							phaseFactorInc += twoPowerLevel;
						}//nn
					}// index
				}//pp
		    	
		    	/*
		    	 * multiply in frequency domain
		    	 */
		    	
		    	for (int nn = 0; nn < N; nn++) {
		    		float aTmp = currentFdVectorReal[nn];
		    		
		    		//real part
		    		currentFdVectorReal[nn] *= codesDFTReal[nn];
		    		currentFdVectorReal[nn] -= codesDFTImag[nn]*currentFdVectorImag[nn]; 
		    		
		    		//imaginary part
		    		currentFdVectorImag[nn] *= codesDFTReal[nn];  
		    		currentFdVectorImag[nn] += codesDFTImag[nn]*aTmp;
		    		currentFdVectorImag[nn] *= -1;  //apply complex conjugation to use forward fft as inverse
		    	}// element-wise multiplication in frequency domain
		    	
		    	/*
		    	 * Radix-2 DIT FFT of sample vector -> works as inverse fft 
		    	 */
		    		//final int stages = nExponent2Final; yet defined
				final float phaseFactorBase = minus2PI/N;
				
				for (int pp = 0; pp < stages; pp++) {
					final int twoPowerLevel = 2<<(pp);
					final int offset = 1<<(pp);
					
					final int phaseFactorInc = (N>>(pp+1));
					for (int index = 0; index < N; index += twoPowerLevel) {
						int kk = 0;
						for (int nn = index; nn < index+offset; nn++) {
							//angle
							final float phaseFactor  = phaseFactorBase*kk;
							
							//phase factor values re, im
							final float phaseFactorReal = PETrigonometry.cos(phaseFactor);
							final float phaseFactorImag = PETrigonometry.sin(phaseFactor);
							
							//multiply with phase factor
							  float b0Real  = phaseFactorReal*currentFdVectorReal[nn+offset];
									b0Real -= phaseFactorImag*currentFdVectorImag[nn+offset];
									
							  float b0Imag  = phaseFactorReal*currentFdVectorImag[nn+offset];
							        b0Imag += phaseFactorImag*currentFdVectorReal[nn+offset];
									
							//assign to memory
							        currentFdVectorReal[nn+offset] = currentFdVectorReal[nn] - b0Real;
							        currentFdVectorImag[nn+offset] = currentFdVectorImag[nn] - b0Imag;
							
							currentFdVectorReal[nn]        = currentFdVectorReal[nn] + b0Real;
							currentFdVectorImag[nn]			= currentFdVectorImag[nn] + b0Imag;
							kk += phaseFactorInc;
						}//nn
					}//index
					//if (twoPowerLevel==N) {
					//	System.out.println("break in startAqc");
					//	break;
					//}
				}//pp
	    		
	    		/*
		    	 *  apply time domain aliasing to rearrange tau indices correctly 
		    	 *  (downsampling in frequency domain, but more accurate in time domain in our case)
		    	 */
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	int maxTauIndexFd = 0;
		    	int ii;
		    	for (ii = 0; ii <nrOfSamples; ii++) { 
		    		//aliasing: add up what should have been wrapped
		    		currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
		    		currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
		    		currentFdVectorReal[ii] /= N;
		    		currentFdVectorImag[ii] /= N;
		    		
		    		float tmp  = currentFdVectorReal[ii]*currentFdVectorReal[ii];
		    		float tmp2 = currentFdVectorImag[ii]*currentFdVectorImag[ii];
		    		
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxTauIndexFd = ii;
		    			currentRFDTimesNMaximum = tmp;
		    			
		    		}// if abs^2 > currentMaximum
		    	}// aliasing and local maximum detection
		    	
	    		maxTauIndices[fd] = maxTauIndexFd;
	    		maxSignalPowersTimesN[fd] = currentRFDTimesNMaximum;
	    	}//fd
	    }//radix 2
	    

		
		
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
	   
	    dopplerShift = maxFdIndex*stepRate + minRate; 
	    acqResult =  maxSignalPowerTimesN > pInTimesNrOfSamples*gamma; //division by N had been omitted in inverse FFT; this way, only needed once
		return acqResult;
	}//startAcquisition
	
	public static void printArraysOctaveFormat(String msg, boolean octave, float[] real, float[] imag) {
		if (msg.length() ==0) return;
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
		final float phaseFactorBase = minus2PI/N;
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
	}//computeCodesDFTIterativelyRadix4
	
	/**
	 * For Testing the inverse DIT approach
	 * Computes the FFT (without scaling 1/N) of the codes vector with an iterative Decimation-in-Time radix-4
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
			
			final int phaseFactorInc = (N>>((pp<<1)+2)); ; //TODO exchange with shifts
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
	   

	/**
	 * Computes the DFT as Radix-2-FFT, iterative, as a Decimation-in-Frequency approach.
	 */
	public void computeCodesDFTIterativelyRadix2DIF() {
		final int stages = nExponent2Final;
		final float phaseFactorExp = minus2PI/N;
		
		for (int pp = 0; pp < stages; pp++) {
			final int twoPowerLevel = 1<<(pp);
			final int halfsize = N>>(pp+1);
			final int size = N>>(pp);
			
			for (int index = 0; index < twoPowerLevel; index++) {
				final int beginKK = index*size;
				final int endKK = beginKK+halfsize;
				int phaseFactorInc = 0;
				for (int nn = beginKK; nn < endKK; nn++) {
					 
					
					final float tmpReal = codesDFTReal[nn];
					final float tmpImag = codesDFTImag[nn];
					
					final float phaseFactor = phaseFactorInc*phaseFactorExp;
					final float phaseFactorReal = PETrigonometry.cos(phaseFactor);
					final float phaseFactorImag = PETrigonometry.sin(phaseFactor);
					
					//System.out.println("twiddlefactor = "+phaseFactorReal+" + "+phaseFactorImag+"i)");
					
					final float a0Real = tmpReal - codesDFTReal[nn+halfsize];
					final float a0Imag = tmpImag - codesDFTImag[nn+halfsize];
					
					float b0Real  = a0Real*phaseFactorReal;
						  b0Real -= a0Imag*phaseFactorImag;
						  
					float b0Imag  = a0Imag*phaseFactorReal;
					      b0Imag += a0Real*phaseFactorImag;
					
					//assignments to memory
					codesDFTReal[nn] = tmpReal + codesDFTReal[nn+halfsize];
					codesDFTImag[nn] = tmpImag + codesDFTImag[nn+halfsize];
					
					codesDFTReal[nn+halfsize] = b0Real;
					codesDFTImag[nn+halfsize] = b0Imag;
					
					phaseFactorInc += twoPowerLevel;
				}//nn
			}// index
		}//pp
	}	//computeCodesDFTIterativelyRadix2DIF

	
	/**
	 * For Testing the inverse DIT approach
	 * Computes the FFT (without scaling 1/N) of the codes vector with an iterative Decimation-in-Time radix-2
	 * approach. 
	 */
	public void computeCodesDFTIterativelyRadix2DIT() {
		final int stages = nExponent2Final;
		final float phaseFactorBase = minus2PI/N;
		
		for (int pp = 0; pp < stages; pp++) {
			final int twoPowerLevel = 2<<(pp);
			final int offset = 1<<(pp);
			
			final int phaseFactorInc = (N>>(pp+1));
			for (int index = 0; index < N; index += twoPowerLevel) {
				int kk = 0;
				for (int nn = index; nn < index+offset; nn++) {
					//angle
					final float phaseFactor  = phaseFactorBase*kk;
					
					//phase factor values re, im
					final float phaseFactorReal = PETrigonometry.cos(phaseFactor);
					final float phaseFactorImag = PETrigonometry.sin(phaseFactor);
					
					//multiply with phase factor
					  float b0Real  = phaseFactorReal*codesDFTReal[nn+offset];
							b0Real -= phaseFactorImag*codesDFTImag[nn+offset];
							
					  float b0Imag  = phaseFactorReal*codesDFTImag[nn+offset];
					        b0Imag += phaseFactorImag*codesDFTReal[nn+offset];
							
					//assign to memory
				    codesDFTReal[nn+offset] = codesDFTReal[nn] - b0Real;
					codesDFTImag[nn+offset] = codesDFTImag[nn] - b0Imag;
					
					codesDFTReal[nn]        = codesDFTReal[nn] + b0Real;
					codesDFTImag[nn]		= codesDFTImag[nn] + b0Imag;
					kk += phaseFactorInc;
				}//nn
			}//index
			//if (twoPowerLevel==N) {
			//	break;
			//}
		}//pp
	}//computeCodesDFTIterativelyRadix2DIT
	
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
