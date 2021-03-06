package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition {
	private static final int sampleThreshold = 50;
	
	//atm just for simulatedAnnealing reasons (composition optimization). 
	public static void main (String[] args) {
		Acquisition acq = new Acquisition(400);
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
		private final static float PI = (float) 3.1415927;
		
		//1 divided by sampling rate
		private final static float oneByFS = 1/400000;
		
		/*
		 * Pre-Computation of frequency shifts
		 *
			// first index: fd; second index: nn
			private  static float[][] precomputedFrequencyShiftsReal;
			private  static float[][] precomputedFrequencyShiftsImag;
			private final static int nrPrecomputedFrequencyShifts = 4096;
		*/
		
		/*
		 * Find valid acquisition
		 */
			private final static float gamma 		= (float)0.015;
		
		/*
		 * samples and codes
		 */
			private final float[] samplesReal;
			private final float[] samplesImag;
		
			private  float[] codesReal;
			private  float[] codesImag;
		
			private float[] codesDFTReal;
			private float[] codesDFTImag;
			
			private float[] codesRefReal;
			private float[] codesRefImag;
			
		
		
		/*
		 *  numbers for samples, FFT
		 */
			private final int nrOfSamples;
			private int nExponent2		= 1;
			private final int N;
			
		
		//utlitiy
			private final float oneDivN;
			private final float oneDivNrOfSamples;
			
			private final boolean fastConvolution;
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
	public Acquisition(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
	    this.oneDivNrOfSamples = 1/nrOfSamples;
	    
	    
	    
	    
		if (nrOfSamples > sampleThreshold) { // use fast convolution via ffts
		    fastConvolution = true;
		    
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
		    oneDivN = (float)1/N;
		    if ((nExponent2 & 0x1) == 1) {
		    	radix4 = false;
		    }
		    else {
		    	radix4 = true;
		    }
		    //
		    
		    //TODO comment out for final hand in
		     System.out.println("nrOfSamples:"+nrOfSamples);
		     System.out.println("N:"+N);
		     System.out.println("exponent:"+nExponent2);
		     System.out.println("radix4 "+ (radix4 ? "possible" : "impossible"));
		    
		    //enteredCodes = N-1;
		    enteredCodes = 0;
		    
		    /*
		     * Note: automatic zero padding due to Java's zero-initial values in arrays.
		     */
		    samplesReal = new float[N];
		    samplesImag = new float[N];
		    
		    codesReal = new float[N];
		    codesImag = new float[N];
		    
		    codesDFTReal = codesReal;
		    codesDFTImag = codesImag;
		    
		    codesRefReal = new float[N];
		    codesRefImag = new float[N];
		    
		}//fast convolution
		
		else { // use classical convolution
			fastConvolution = false;
			radix4 = false;
			
			samplesReal = new float[nrOfSamples];
			samplesImag = new float[nrOfSamples];
			
			codesReal = new float[nrOfSamples];
			codesImag = new float[nrOfSamples];
			
			
			codesDFTReal = null;
			codesDFTImag = null;
			
			 codesRefReal = null;//nullnew float[N];
			 codesRefImag = null;//null new float[N];
			
			N = nrOfSamples;
			oneDivN = oneDivNrOfSamples;
			enteredCodes = nrOfSamples-1;
			
		}//classical convolution
	}
	
	/* not needed by now
	private static void precomputeFrequencyShifts() {
		precomputedFrequencyShiftsReal = new float[nrFrequencies][];
		precomputedFrequencyShiftsImag = new float[nrFrequencies][];
			for (int fd = 0; fd < nrFrequencies; fd++) {
				precomputedFrequencyShiftsReal[fd] = new float[nrPrecomputedFrequencyShifts];
				precomputedFrequencyShiftsImag[fd] = new float[nrPrecomputedFrequencyShifts];
			}
			for (int fd = 0; fd < nrFrequencies; fd++) {
				precomputedFrequencyShiftsReal[fd][0] = 1;
				precomputedFrequencyShiftsImag[fd][0] = 0;
				for (int nn = 1; nn < nrPrecomputedFrequencyShifts; nn++) {
					float phase1 = PI*oneByFS*-2*nn;
					float phase2 = (fd*stepRate + minRate);
									
					float phase = phase1*phase2;
					precomputedFrequencyShiftsReal[fd][nn] = PETrigonometry.cos(phase);
					precomputedFrequencyShiftsImag[fd][nn] = PETrigonometry.sin(phase);
				}// nn loop
			}// fd loop
	}//precomputeFrequencyShifts
	*/
	
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
		//int index = indexPointersBitReversed[enteredCodes];
		// normal codes
		codesDFTReal[enteredCodes] = real;
		codesDFTImag[enteredCodes] = imag;
		codesRefReal[enteredCodes] = real;
		codesRefImag[enteredCodes] = imag;
		//enteredCodes--;
		enteredCodes++;
	} //enterCode
	
	public boolean startAcquisition(){
		int maxFdIndex = 0;
	    int maxTauIndex = 0;
	    float maxSignalPowerTimesN = 0;
	    int[] maxTauIndices = new int[nrOfSamples];
	    float[] maxSignalPowersTimesN = new float[nrFrequencies];
	    boolean acqResult = false;
	    
	    /*
	     * nrOfSamples is such that a fast convolution, as applied in that case, may be replaced by a radix-4 version, which yields significant speedup
	     * over standard cyclic convolution and over radix 2 (25% less comlex multiplications)
	     */
	    if (radix4) {
	    	float[] currentFdVectorReal = new float[N];
		    float[] currentFdVectorImag = new float[N];
		    final float nFactorRadix4 = 1/(N);
		    final float nFactorRadix4Squared = nFactorRadix4*nFactorRadix4;
		    /*
		     * compute DFT of code vector
		     */
		    
		   // compareCodes();
		   computeCodesDFTIterativelyRadix4DIF();
		   
		   //computeCodesDFTIterativelyRadix4Inv();
		   //compareCodes();
		   
		   
		    //System.exit(0);
		    /*
		     * Loop: Do the following for each frequency 
		     */
		   final int stages = nExponent2>>1;
		   final float phaseFactorBase = -2*PI*nFactorRadix4;
		   final float phaseFactorBaseInverse = -phaseFactorBase;
		   for (int fd = 0; fd < nrFrequencies; fd++) {
			   
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
		    	
			   
		    	/*
		    	 * Radix 4 FFT of Samples for current Frequency
		    	 */
		    	
				
				
			   
		    	
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
		    	 * Inverse Radix 4 FFT for current frequency
		    	 */
		    	
		    	
				
			   
			   
			   
			   /*
		    	 *  apply time domain aliasing to rearrange tau indices correctly 
		    	 *  (downsampling in frequency domain, but more accurate in time domain in our case)
		    	 *       apply scaling by 1/N after inverse FFT here, also
		    	 */
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	int ii = 0;
		    	for (; ii <nrOfSamples-1; ) { //TODO check how to implement maximum check for last step 
		    		currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
		    		currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
		    		float tmp = currentFdVectorReal[ii]*currentFdVectorReal[ii];
		    		float tmp2= currentFdVectorImag[ii]*currentFdVectorImag[ii];
		    		tmp  *= nFactorRadix4Squared;
		    		tmp2 *= nFactorRadix4Squared;
		    		
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxSignalPowersTimesN[fd] = tmp;
		    			maxTauIndices[fd] = ii;
		    			currentRFDTimesNMaximum = tmp;
		    			
		    		}// if abs^2 > currentMaximum
		    		ii++;
		    	}// aliasing and local maximum detection
		    	//last one (unaliased)
	    		float tmp = currentFdVectorReal[ii]*currentFdVectorReal[ii];
	    		float tmp2= currentFdVectorImag[ii]*currentFdVectorImag[ii];
	    		tmp  *= nFactorRadix4Squared;
	    		tmp2 *= nFactorRadix4Squared;
	    		
	    		if (tmp>currentRFDTimesNMaximum-tmp2) {
	    			tmp = tmp+tmp2;
	    			maxSignalPowersTimesN[fd] = tmp;
	    			maxTauIndices[fd] = ii;
	    			currentRFDTimesNMaximum = tmp;
	    			
	    		}// if abs^2 > currentMaximum
		   }//fd
	    } //radix 4
	    
	    /*
	     * nrOfSamples is so high that fast convolution brings advantage over classical
	     */
	    else if (fastConvolution) {

		    float[] currentFdVectorReal = new float[N];
		    float[] currentFdVectorImag = new float[N];
		    /*
		     * compute DFT of code vector
		     */
		   computeCodesDFTIteratively();
		    
		    /*
		     * Loop: Do the following for each frequency 
		     */
		    for (int fd = 0; fd < nrFrequencies; fd++) {
		    	
		    	//final int frequency = minRate+fd*stepRate;
		    	final float fdEFactorPhase = (float)-0.01570796326*fd; // = (-2*PI*stepRate*1/fs) * fd
		    	final float fdPrecomputeFactorPhase = fdEFactorPhase + fixPrecomputeEFactorPhase;
		    	
		    	            //in this iteration, work on the precomputed factors for the samples of frequency[fd]
		    	            	//currentFdVectorReal = precomputedFrequencyShiftsReal[fd];
		    	            	//currentFdVectorImag = precomputedFrequencyShiftsImag[fd];
		    	
		    		
		    		//currentFdVectorReal[0] = 1;
		    		//currentFdVectorImag[0] = 0; // can be omitted due to java zero initialization
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
		    	
		    	            /* not necessary any more; keep for understanding what is actually happening (zero padding for fast convolution)
		    	            //artificial padding, because its faster than multiplication with zero
		    	            for (int nno = nrOfSamples; nno < N; nno++) {
		    	            	//apply zero padding afterwards 
		    	            	currentFdVectorReal[nno] = 0;
		    	            	currentFdVectorImag[nno] = 0;
		    	            }
	        	            */
		    	
		    	
		    	/*
		    	 *  computation of fft of multiplied samples for current fd frequency
		    	 */
		    	
		    	int p = 1;
		    	for (int ii = 2; ii <=N; ii*=2) { //fft outermost loop: recursion level
		    		
		    		final float phaseConst = -2*(1/ii)*PI; 
		    		//for kk=0, skip long computation of phase
		    		final int nByII = N>>>p;
		    		p++;
		    		
		    		
		    		//then, the original loop
		    		for (int kk = 1; kk < ii/2; kk++) {
		    			//compute phase
		    			final float phase = phaseConst*kk;
		    			final float phaseReal = PETrigonometry.cos(phase);
		    			final float phaseImag = PETrigonometry.sin(phase);
		    			
		    			for (int jj = 0; jj < nByII; jj++) {
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
		    	p = 1;
		    	for (int ii = 2; ii <=N; ii*=2) { //fft outermost loop: recursion level
		    		
		    		final float phaseConst = -2*(1/ii)*PI; 
		    		//for kk=0, skip long computation of phase
		    		final int nByII = N>>>p;
		    		p++;
		    		
		    		
		    		//then, the original loop
		    		for (int kk = 1; kk < ii/2; kk++) {
		    			//compute phase
		    			final float phase = phaseConst*kk;
		    			final float phaseReal = PETrigonometry.cos(phase);
		    			final float phaseImag = PETrigonometry.sin(phase);
		    			
		    			for (int jj = 0; jj < nByII; jj++) {
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
		    	 *  apply time domain aliasing to rearrange tau indices correctly 
		    	 *  (downsampling in frequency domain, but more accurate in time domain in our case)
		    	 */
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	for (int ii =0; ii <nrOfSamples-1; ii++) { //TODO check how to implement maximum check for last step 
		    		currentFdVectorReal[ii] += currentFdVectorReal[ii+nrOfSamples];
		    		currentFdVectorImag[ii] += currentFdVectorImag[ii+nrOfSamples];
		    		float tmp = currentFdVectorReal[ii]*currentFdVectorReal[ii];
		    		float tmp2= currentFdVectorImag[ii]*currentFdVectorImag[ii];
		    		
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxSignalPowersTimesN[fd] = tmp;
		    			maxTauIndices[fd] = ii;
		    			currentRFDTimesNMaximum = tmp;
		    			
		    		}// if abs^2 > currentMaximum
		    		
		    	}// aliasing and local maximum detection
		    	
		    } //fd
		    
   
		}//fast convolution
		
		/*
		 * classical, circular convolution
		 */
		else { 
			
			float[] currentFdVectorReal = new float[nrOfSamples];
		    float[] currentFdVectorImag = new float[nrOfSamples];
		    
		    float[] circularResultReal = new float[nrOfSamples];
		    float[] circularResultImag = new float[nrOfSamples];

		    /*
		     * Loop: Do the following for each frequency 
		     */
		    for (int fd = 0; fd < nrFrequencies; fd++) {
		    	
		    	//final int frequency = minRate+fd*stepRate;
		    	final float fdEFactorPhase = (float)-0.01570796326*fd; // = (-2*PI*stepRate*1/fs) * fd
		    	final float fdPrecomputeFactorPhase = fdEFactorPhase + fixPrecomputeEFactorPhase;
		    			    		
		    		//currentFdVectorReal[0] = 1;
		    		//currentFdVectorImag[0] = 0; // can be omitted due to java zero initialization
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
		    	
		    	/*
		    	 * circular convolution
		    	 */
		    	for (int nn = 0; nn < nrOfSamples; nn++) {
		    		

		    		// index is the evaluation of (n-i)mod nrOfSamples
		    		//   for each loop iteration of the outer loop (nn = const)
		    		//   index starts a nn, then goes down to zero in the first loop,
		    		//   then goes up again; in total, index has nrOfSamples values in one
		    		//   iteration of the outer loop
		    		int index1 = nn;
		    		int index2 = -nn;
		    		
		    		float accReal = 0;
		    		float accImag = 0;
		    		float accReal2 = 0;
		    		float accImag2 = 0;
		    		for (int ii = 0; ii<nrOfSamples; ii++) {
		    			//multiplication
		    			float factorAReal = currentFdVectorReal[ii];
		    		    float factorAImag = currentFdVectorImag[ii];
		    		    int index = ii < nn ? index1 : index2;
		    		    float factorBReal = codesReal[index];
		    		    float factorBImag = codesImag[index];
		    		    
		    		    accReal += factorAReal*factorBReal;
		    		    accReal -= factorAImag*factorBImag;
		    		    accImag += factorAReal*factorBImag;
		    		    accImag += factorBReal*factorAImag;
		    		    
		    		    // adapt index;
		    			index1--;
		    			index2++;
		    		}// ii
		    		

		    		
		    		
		    		circularResultReal[nn] = accReal;
		    		circularResultImag[nn] = accImag;
		    	}// nn
		    	
		    	float currentRFDTimesNMaximum = 0;
		    	for (int ii =0; ii <nrOfSamples; ii++) {
		    		float tmp = circularResultReal[ii]*circularResultReal[ii];
		    		float tmp2= circularResultImag[ii]*circularResultImag[ii];
		    		if (tmp>currentRFDTimesNMaximum-tmp2) {
		    			tmp = tmp+tmp2;
		    			maxSignalPowersTimesN[fd] = tmp;
		    			maxTauIndices[fd] = ii;
		    			currentRFDTimesNMaximum = tmp;			
		    		}// if abs^2 > currentMaximum
		    	}// local maximum detection
		    }//fd
			
		}//classical convolution
		
		
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
	    acqResult =  maxSignalPowerTimesN > pInTimesNrOfSamples*gamma;
		return acqResult;
	}//startAcquisition
	
	public void computeCodesDFTIterativelyRadix4DIF() {
		final int stages = nExponent2>>1; //log(2)/2
		final float phaseFactorBase = -2*PI/N;
		int M = N>>2; // N/4
		
		//outermost loop: recursion level
		for (int pp = 0; pp<stages; pp++) {
			
			int fourPowerLevel = 1 << (pp<<1);
			
			
			float phaseFactorP1Base = 2*fourPowerLevel*phaseFactorBase;
			float phaseFactorP2Base = fourPowerLevel*phaseFactorBase;
			float phaseFactorP3Base = 3*fourPowerLevel*phaseFactorBase;
			
			// middle loop: index of block
			int inc = N>>(pp<<1);
			for (int index = 0; index < N; index+=inc) {
				
				//innermost loop: perform on block
				for (int nn = 0; nn < M; nn++) {
					
					int indexM = index+M; // testing if faster if in outer loop
					int index2M = index+2*M; //see line above
					int index3M = index+3*M; //see line above
					
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
						resultImag2 = b0Imag-b1Imag-b2Imag+b3Real;
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
						resultImag3 = b0Imag+b1Imag-b2Imag-b3Real;
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


	public void computeCodesDFTIterativelyRadix4DITInverse() {
		final int stages = nExponent2>>1; // log2(N)/2;
		
		final float phaseFactorBase = 2*PI/N; //no minus because of inversion
		
		for (int pp = 0; pp<stages; pp++) {
			int fourPowerLevel = 4<<(pp<<1);
			int offset = 1<<(pp<<1);
			int offset2 = offset<<1;
			int offset3 = 3*offset;
			int phaseFactorInc = N/fourPowerLevel; //TODO excahnge with shifts
			for (int index = 0; index < N; index+=fourPowerLevel) {
				int kk = 0;
				for (int nn = index; nn < index+offset; nn++) {
					
					//phase factors
					float phaseFactorKK = phaseFactorBase*kk;
					
					float phaseFactor1Real = PETrigonometry.cos(phaseFactorKK);
				    float phaseFactor1Imag = PETrigonometry.sin(phaseFactorKK);
				    
				    float phaseFactor2Real = PETrigonometry.cos(2*phaseFactorKK);
				    float phaseFactor2Imag = PETrigonometry.sin(2*phaseFactorKK);
				    
				    float phaseFactor3Real = PETrigonometry.cos(3*phaseFactorKK);
				    float phaseFactor3Imag = PETrigonometry.sin(3*phaseFactorKK);
				    		
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
						//complex a0 = b0*phaseFactor0 (1);
						//complex a1 = b2*phaseFactor1 
					    //complex a2 = b1*phaseFactor2
					    //complex a3 = b3*phaseFactor3;#
					
					float a0Real, a0Imag;
					float a1Real, a1Imag;
					float a2Real, a2Imag;
					float a3Real, a3Imag;
					
					//a0
					a0Real = b0Real;
					a0Imag = b0Imag;
					
					//a1
					a1Real = b2Real*phaseFactor1Real - b2Imag*phaseFactor1Imag;
					a1Imag = b2Real*phaseFactor1Imag + b2Imag*phaseFactor1Real;
					
					//a2
					a2Real = b1Real*phaseFactor2Real - b1Imag*phaseFactor2Imag;
					a2Imag = b1Real*phaseFactor2Imag + b1Imag*phaseFactor2Real;
					
					//a3
					a3Real = b3Real*phaseFactor3Real - b3Imag*phaseFactor3Imag;
					a3Imag = b3Real*phaseFactor3Imag + b3Imag*phaseFactor3Real;
					
					//assign to memory
					codesDFTReal[nn] = a0Real + a2Real + a1Real + a3Real;
					codesDFTImag[nn] = a0Imag + a2Imag + a1Imag + a3Imag;
					
					codesDFTReal[nn+offset] = a0Real - a2Real;
					codesDFTImag[nn+offset] = 0;
					
					codesDFTReal[nn+offset2] = 0;
					codesDFTImag[nn+offset2] = 0;
					
					codesDFTReal[nn+offset3] = 0;
					codesDFTImag[nn+offset3] = 0;
				    		
					kk +=phaseFactorInc;
				}//nn loop
			}//index
		}//pp = recursion level
		
	}//computeCodesDFTIterativelyRadix4DIT
	   
	public void compareCodes() {
		float byN = 1/N;
		
		//aliasing ? 
		for (int nn=0; nn < nrOfSamples-1; nn++) {
			codesDFTReal[nn] +=codesDFTReal[nn+nrOfSamples];
			codesDFTImag[nn] +=codesDFTImag[nn+nrOfSamples];
		}
		for (int nn = 0; nn < N; nn++) {
			double realDiff = codesDFTReal[nn]/N - codesRefReal[nn];
			double imagDiff = codesDFTImag[nn]/N - codesRefImag[nn];
			String line = nn+ ":\tref: ("+codesRefReal[nn]+","+codesRefImag[nn]+")  \t computed:("+codesDFTReal[nn]+","+codesDFTImag[nn]+")\t diff: ("+realDiff+","+imagDiff+")";
			System.out.println(line);
		}
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
 * Getters
 */

public int getDopplerverschiebung(){
	return dopplerShift;
}//getDopplerverschiebung

public int getCodeVerschiebung(){
	return codeShift;
} // getCodeVerschiebung

}
