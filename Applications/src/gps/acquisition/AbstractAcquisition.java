package gps.acquisition;

public abstract class AbstractAcquisition {
	
	protected final static int samplingRate = 400000; //400kHz
	protected final static int stepRate 	= 1000; //1 kHz
	protected final static int maxRate	 	= 5000; //5kHz
	protected final static int minRate 		= -5000; //-5kHz
	
	protected final static int nrFrequencies = Math.floorDiv((maxRate-minRate), stepRate)+1;
	
	protected final static int[] allFrequencies = computeAllFrequencies();
	
	protected final static float[][][] precomputedFrequencyShifts = computeFrequencyShifts();
	
	protected final static float gamma 		= (float)0.015;
	
	protected final static int nDivisor		= 16;
	
	protected final static int maxFactorPreCompute = 256;
	
	protected final int nrOfSamples;
	protected final float oneDivN;
	
	protected boolean convo;

	public AbstractAcquisition(int nrOfSamples){
		this.nrOfSamples = nrOfSamples;
		this.oneDivN = (float) 1/nrOfSamples;
		convo = false;
	}
	
	public void setConvo(boolean convo) {
		this.convo = convo;
	}
	
	/**
	 * Pre-computes all existing frequencies
	 * @return an int[] with all the frequencies
	 */
	private static int[] computeAllFrequencies() {
		int[] f = new int[nrFrequencies];
		for (int fd = 0; fd < nrFrequencies; fd++) {
			f[fd] = minRate+(fd*stepRate);
		}
		return f;
	}
	
	private static float[][][] computeFrequencyShifts(){
		float[][][] shifts = new float[nrFrequencies][][];
		for (int fd = 0; fd< nrFrequencies; fd++) {
			shifts[fd] = new float[nDivisor*maxFactorPreCompute][];
			for (int nn = 0; nn < nDivisor*maxFactorPreCompute; nn++) {
				shifts[fd][nn] = getCartesianFormFromPolar(1, (float)-2*(float)Math.PI*nn*allFrequencies[fd]/samplingRate);
			} // n
		} // m
		return shifts;
	}

	public abstract void enterSample(float real, float imag);
	
	public abstract void enterCode(float real, float imag);
	
	public abstract boolean startAcquisition();
	
	public abstract int getDopplerverschiebung();
	
	public abstract int getCodeVerschiebung();
	
	public static float[] getCartesianFormFromPolar(float abs, float phase) {
		float[] res = new float[2];
		res[0] = (float)Math.cos(phase)*abs;
		res[1] = (float)Math.sin(phase)*abs;
		return res;
	}
	
	public static float[] cartesianComplexMult(float[] a, float[] b) {
		float[] result = new float[2];
		result[0] = cartesianComplexMultReal(a, b);
		result[1] = cartesianComplexMultImag(a, b);
		return result;
	}
	
	public static float cartesianComplexMultReal(float[] a, float[] b) {
		float result = a[0]*b[0]-a[1]*b[1];
		return result;
	}
	
	public static float cartesianComplexMultImag(float[] a, float[] b) {
		float result = a[0]*b[1]+a[1]*b[0];
		return result;
	}
	
	public static float[] polarComplexMult(float[] a, float[] b) {
		float[] result = new float[2];
		result[0] = polarComplexMultAbs(a[0], b[0]);
		result[1] = polarComplexMultPhase(a[1], b[1]);
		return result;
	}
	
	public static float polarComplexMultAbs(float a, float b) {
		return a*b;
	}
	
	public static float polarComplexMultPhase(float a, float b) {
		return a+b;
	}
	
	public static float polarComplexMultAbs(float[] a, float[] b) {
		return a[0]*b[0];
	}
	
	public static float polarComplexMultPhase(float[] a, float[] b) {
		return a[1]+b[1];
	}
	
	public static void polarComplexMultConst(float[] complex, float constant) {
		complex[0] *= constant;
	}
	
	public static void cartesianComplexMultConst(float[] complex, float constant) {
		complex[0] *= constant;
		complex[1] *= constant;
	}
	
	
	public static void cartesianComplexAccumulate(float[] sourceDest, float[] toAdd) {
		sourceDest[0] += toAdd[0];
		sourceDest[1] += toAdd[1];
	}
	
	public static float[][] complexCyclicConvolution(float[][] f, float[][] g){
		int N = f.length;
		float[][] convResult = new float[N][];
		for (int nn = 0; nn < N; nn++) {
			float[] currentNValue = new float[] {0, 0};
			for (int kk = 0; kk <=nn ; kk++) {
				cartesianComplexAccumulate(	currentNValue, 
											cartesianComplexMult(f[kk], g[(nn-kk)%N]));	//original: nn-kk
			}
			convResult[nn] = currentNValue;
		}
		return convResult;
	}
	
	public static float[][] complexDFT(float[][] vector){
		final int N = vector.length;
		final float oneByN = (float)1/N;
		float[][] dft = new float[N][];
		for (int kk = 0; kk < N; kk++) {
			dft[kk] = new float[] {0, 0};
			for (int nn = 0; nn < N; nn++) {
				float[] toAdd = cartesianComplexMult(	vector[nn], 
														getCartesianFormFromPolar(1, -2*(float)Math.PI*nn*kk*oneByN));
				cartesianComplexAccumulate(dft[kk], toAdd);
			}
		}
		return dft;
	}
	
	public static float[][] complexVectorMult(float[][] a, float[][] b){
		final int N = a.length;
		float[][] prod = new float[N][];
		for (int nn = 0; nn < N; nn++) {
			prod[nn] = cartesianComplexMult(a[nn], b[nn]);
		}
		return prod;
	}
	
	public static float[][] complexIDFT(float[][] vector){
		final int N = vector.length;
		final float oneByN = (float)1/N;
		float[][] idft = new float[N][];
		for (int nn = 0; nn < N; nn++) {
			
			float[] acc = new float[] {0, 0};
			for (int kk = 0; kk < N; kk++) {
				float[] toAdd = cartesianComplexMult(	vector[kk],
														getCartesianFormFromPolar(1, 2*(float)Math.PI*kk*nn*oneByN));
				cartesianComplexAccumulate(	acc, toAdd);
			}
			
			idft[nn] = acc;
			cartesianComplexMultConst(acc, oneByN);
		}
		return idft;
	}
	
	public static float[][] transponeVector(float[][] vector){
		final int N = vector.length;
		float[][] trans = new float[N][];
		for (int nn = 0; nn < N; nn++) {
			trans[N-1-nn] = vector[nn];
		}
		return trans;
	}
	
	public static float[][] cartesianComplexConjugateVector(float[][] vector){
		float[][] conj = new float[vector.length][];
		for (int ii = 0; ii < vector.length; ii++) {
			conj[ii] = new float[] {vector[ii][0], -vector[ii][1]};
		}
		return conj;
	}
}
