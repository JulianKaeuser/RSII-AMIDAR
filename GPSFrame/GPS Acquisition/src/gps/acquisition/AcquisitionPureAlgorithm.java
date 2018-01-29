package gps.acquisition;

public class AcquisitionPureAlgorithm extends AbstractAcquisition {
	
	//sample related
	private int nrOfSamplesEntered = 0;
	private float[][] samples;
	
	//code related
	private int nrOfCodesEntered = 0;
	private float[][] codes;
	
	private float[][] complexDFTCodesConj;
	
	
	// average power
	private float pInTimesN;
	
	//results
	private int codeShift;
	private int dopplerShift;
	private boolean acquisitionResult;

	/**
	 * Constructor, taking the nr of samples as argument
	 * @param nrOfSamples
	 */
	public AcquisitionPureAlgorithm(int nrOfSamples) {
		super(nrOfSamples);
		samples = new float[nrOfSamples][];
		codes   = new float[nrOfSamples][];

		pInTimesN = 0;
		codeShift = 0;
		dopplerShift = 0;
	}

	public void enterSample(float real, float imag) {
		if(nrOfSamplesEntered<nrOfSamples) {
			samples[nrOfSamplesEntered] = new float[] {real, imag};
			pInTimesN += (real*real)+(imag*imag);	
			nrOfSamplesEntered++;
		}	
	}

	public void enterCode(float real, float imag) {
		if (nrOfCodesEntered<nrOfSamples) {
			codes[nrOfCodesEntered] = new float[] {real, imag};
			nrOfCodesEntered++;
		}	
	}

	/* #####################################################
	 * ########### Algorithm ###############################
	 * #####################################################
	 */
	
	/**
	 * Starts the algorithm to find the acquisition
	 */
	public boolean startAcquisition() {
		float[][] codesTrans = transponeVector(codes);
		float[][] codesConj = cartesianComplexConjugateVector(codes);
		float[][] codesConjTrans = transponeVector(codesConj);
		//float[][] codesConjTrans = codesConj;
		complexDFTCodesConj = complexDFT(codesTrans); //working: codesTrans, codeConjTrans
		
		
		int maxFdIndex = 0;
		int maxTauIndex = 0;
		float maxSignalPowerTimesN = 0;
		float[][] currentRFDTimesN;
		
		/*
		 *find Xfd values
		 */
		
		//current Xfd-line is reused every time; no "real" matrix operations necessary, therefore, line by line
		float[][] currentXMatrixLine = new float[nrOfSamples][];
		
		// for each frequency, do:
		for (int fd = 0; fd < nrFrequencies; fd++) {
			
			// multiply input samples nn with pre-computed e^(-j2pi*fd*nn/fsampling)
			for (int nn = 0; nn < nrOfSamples; nn++) {
				currentXMatrixLine[nn] = cartesianComplexMult(samples[nn], precomputedFrequencyShifts[fd][nn]);
			} //nn = running index of sample
			
			// compute IDFT[DFT[currentXMatrixLine]*[DFT[Codes]]
			
			if (convo) {
				currentRFDTimesN = complexCyclicConvolution(currentXMatrixLine, codesConjTrans);
			}
			else {currentRFDTimesN = complexIDFT(	complexVectorMult(
												complexDFT(currentXMatrixLine),
												complexDFTCodesConj));
			}
			
			float currentRFDTimesNMaximum = 0;
			int   currentMaxTauIndex = 0;
			for (int RFDIndex = 0; RFDIndex < currentRFDTimesN.length; RFDIndex++) {
				float tmp = (currentRFDTimesN[RFDIndex][0]*currentRFDTimesN[RFDIndex][0]
						+ currentRFDTimesN[RFDIndex][1]*currentRFDTimesN[RFDIndex][1]);
				if (tmp>currentRFDTimesNMaximum) {
					currentRFDTimesNMaximum = tmp;
					currentMaxTauIndex = RFDIndex;
				}// if abs^2 > currentMaximum
			} //RFDIndex - finding the maximum
			if (currentRFDTimesNMaximum > maxSignalPowerTimesN) {
				maxSignalPowerTimesN = currentRFDTimesNMaximum;
				maxFdIndex = fd;
				maxTauIndex = currentMaxTauIndex;
			}
		}// fd
		
		// define data to return
		codeShift = maxTauIndex+1;
		dopplerShift = allFrequencies[maxFdIndex];
		acquisitionResult = maxSignalPowerTimesN > pInTimesN*gamma ? true : false;
		return acquisitionResult;
	}
	
	/* ####################################################
	 * ######## Getters ###################################
	 * ####################################################
	 */
	
	/*
	 * returns the doppler shift
	 * (non-Javadoc)
	 * @see gps.acquisition.AbstractAcquisition#getDopplerverschiebung()
	 */
	public int getDopplerverschiebung() {
		return dopplerShift;
	}

	public int getCodeVerschiebung() {
		return codeShift;
	}
	
	

}
