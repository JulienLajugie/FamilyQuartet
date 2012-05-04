package dataStructures;

import java.util.List;


/**
 * This class computes the dominant state for each variant using either a moving average algorithm
 * or a binning algorithm
 * @author Julien Lajugie
 */
public class VariantListAnalyzer {

	private static final int HALF_WINDOW_SIZE = 500000;			// half size of the window used for the moving average or for the bins
	private static final boolean COMPUTE_MOVING_WINDOW = false;	// if true compute the inheritance state using a moving window (uses bins otherwise)
	
	private final List<Variant>		variantList;				// list of variants for a family quartet
	private final double[] 			identicals;					// list of identical variants
	private final double[] 			identicalAvgs;				// moving average of the identical variants
	private final double[] 			maternalIdenticals;			// list of haploidentical maternal variants
	private final double[] 			maternalIdenticalAvgs;		// moving average of the haploidentical maternal variants
	private final double[] 			paternalIdenticals;			// list of haploidentical paternal variants
	private final double[]			paternalIdenticalAvgs;		// moving average of the haploidentical paternal variants
	private final double[] 			nonIdenticals;				// list of non identical variants
	private final double[] 			nonIdenticalAvgs;			// moving average of the non identical variants
	private final double[]			nonInfomatives;				// list non informative variants
	private final double[]			nonInfomativeAvgs;			// moving average on non informative variants
	private final double[]			nonMendelians;				// list of non mendelian variants
	private final double[]			nonMendelianAvgs;			// moving average on mendelian variants
	private final QuartetInheritanceState[]dominantInheritanceStates;	// dominant inheritance state of the variants


	/**
	 * Creates an instance of {@link VariantListAnalyzer}
	 * @param variantList 
	 */
	public VariantListAnalyzer(List<Variant> variantList) {
		this.variantList = variantList;
		this.identicals = computeStateScores(QuartetInheritanceState.IDENTICAL);
		this.maternalIdenticals = computeStateScores(QuartetInheritanceState.MATERNAL);
		this.paternalIdenticals = computeStateScores(QuartetInheritanceState.PATERNAL);
		this.nonIdenticals = computeStateScores(QuartetInheritanceState.NON_IDENTICAL);
		this.nonInfomatives = computeStateScores(QuartetInheritanceState.NOT_INFORMATIVE);
		this.nonMendelians = computeStateScores(QuartetInheritanceState.MIE);

		this.identicalAvgs = computeAvgStateScores(identicals);
		this.maternalIdenticalAvgs = computeAvgStateScores(maternalIdenticals);
		this.paternalIdenticalAvgs = computeAvgStateScores(paternalIdenticals);
		this.nonIdenticalAvgs = computeAvgStateScores(nonIdenticals);
		this.nonInfomativeAvgs = computeAvgStateScores(nonInfomatives);
		this.nonMendelianAvgs = computeAvgStateScores(nonMendelians);

		this.dominantInheritanceStates = computeDominantInheritanceStates();
	}
	
	
	/**
	 * Computes the scores of the specified state
	 * @param quartetInheritanceState an inheritance state
	 * @return an array containing a score for each variant for the specified inheritance state
	 */
	private double[] computeStateScores(QuartetInheritanceState quartetInheritanceState) {
		double[] resultState = new double[variantList.size()];
		for (int i = 0; i < variantList.size(); i++) {
			Variant currentVariant = variantList.get(i);
			if (currentVariant.getInheritanceStates()[0] == quartetInheritanceState) {
				double score = 1 / (double) currentVariant.getInheritanceStates().length;
				resultState[i] = score;
			} else if (currentVariant.getInheritanceStates().length > 1) {
				if (currentVariant.getInheritanceStates()[1] == quartetInheritanceState) {
					resultState[i] = 0.5;
				}
			}			
		}
		return resultState;
	}


	/**
	 * Computes the average scores of the specified state.  The average can be computed on bins or on a moving window
	 * @param stateScores array of double containing the scores on which the average algorithm needs to be applied
	 * @return an array containing an averaged score for each variant
	 */
	private double[] computeAvgStateScores(double[] stateScores) {
		if (COMPUTE_MOVING_WINDOW) {
			return computeMovingWindowAvg(stateScores);
		} else {
			return computeAvgStatesBinned(stateScores);
		}
	}


	/**
	 * Computes a moving window average on the specified array of double
	 * @param stateScores array of double
	 * @return an array containing the averaged scores
	 */
	private double[] computeMovingWindowAvg(double[] stateScores) {
		double[] resultAvg = new double[stateScores.length];
		int startIndex = 0;
		
		for (int i = 0; i < stateScores.length; i++) {
			Variant currentVariant = variantList.get(i);
			while ((startIndex < variantList.size()) &&
					(variantList.get(startIndex).getPosition() < currentVariant.getPosition() - HALF_WINDOW_SIZE)) {
				startIndex++;
			}
			int j = startIndex;
			double sum = 0;
			int count = 0;
			while ((j < variantList.size()) &&
					(variantList.get(j).getPosition() <= currentVariant.getPosition() + HALF_WINDOW_SIZE)) {
				sum += stateScores[j];
				count++;
				j++;
			}
			resultAvg[i] = sum / (double) count;			
		}		
		return resultAvg;
	}


	/**
	 * Computes a binned average on the specified array of double
	 * @param stateScores array of double
	 * @return an array containing the averaged scores
	 */
	private double[] computeAvgStatesBinned(double[] stateScores) {
		double[] resultAvg = new double[stateScores.length];
		int currentWindowStop = HALF_WINDOW_SIZE * 2;
		int indexWindowStart = 0; 
		int currentIndex = 0;
		while (currentIndex < variantList.size()) {
			indexWindowStart = currentIndex;
			double sum = 0;
			int count = 0;
			while ((currentIndex < variantList.size()) && 
					(variantList.get(currentIndex).getPosition() < currentWindowStop)) {
				sum += stateScores[currentIndex];
				count++;
				currentIndex++;
			}
			double currentScore = sum / (double) count; 
			for (int i = indexWindowStart; i < currentIndex; i++) {
				resultAvg[i] = currentScore;
			}
			currentWindowStop += HALF_WINDOW_SIZE * 2;
		}
		return resultAvg;
	}
	
	
	/**
	 * @return an array with the dominant inheritance state of each variant.  
	 */
	private QuartetInheritanceState[] computeDominantInheritanceStates() {
		QuartetInheritanceState[] dominantInheritanceStates = new QuartetInheritanceState[variantList.size()];
		for (int i = 0; i < variantList.size(); i++) {
			QuartetInheritanceState dominantState = QuartetInheritanceState.IDENTICAL;
			double maxScore = identicalAvgs[i];
			if (maternalIdenticalAvgs[i] > maxScore) {
				maxScore = maternalIdenticalAvgs[i];
				dominantState = QuartetInheritanceState.MATERNAL;
			}
			if (paternalIdenticalAvgs[i] > maxScore) {
				maxScore = paternalIdenticalAvgs[i];
				dominantState = QuartetInheritanceState.PATERNAL;
			}
			if (nonIdenticalAvgs[i] > maxScore) {
				maxScore = nonIdenticalAvgs[i];
				dominantState = QuartetInheritanceState.NON_IDENTICAL;
			}
			if (nonInfomativeAvgs[i] > maxScore) {
				maxScore = nonInfomativeAvgs[i];
				dominantState = QuartetInheritanceState.NOT_INFORMATIVE;
			}
			if (nonMendelianAvgs[i] > maxScore) {
				dominantState = QuartetInheritanceState.MIE;
			}
			dominantInheritanceStates[i] = dominantState;			
		}		
		return dominantInheritanceStates;
	}


	/**
	 * @return the variantList
	 */
	public final List<Variant> getVariantList() {
		return variantList;
	}


	/**
	 * @return the identicals
	 */
	public final double[] getIdenticals() {
		return identicals;
	}


	/**
	 * @return the identicalAvgs
	 */
	public final double[] getIdenticalAvgs() {
		return identicalAvgs;
	}


	/**
	 * @return the maternalIdenticals
	 */
	public final double[] getMaternalIdenticals() {
		return maternalIdenticals;
	}


	/**
	 * @return the maternalIdenticalAvgs
	 */
	public final double[] getMaternalIdenticalAvgs() {
		return maternalIdenticalAvgs;
	}


	/**
	 * @return the paternalIdenticals
	 */
	public final double[] getPaternalIdenticals() {
		return paternalIdenticals;
	}


	/**
	 * @return the paternalIdenticalAvgs
	 */
	public final double[] getPaternalIdenticalAvgs() {
		return paternalIdenticalAvgs;
	}


	/**
	 * @return the nonIdenticals
	 */
	public final double[] getNonIdenticals() {
		return nonIdenticals;
	}


	/**
	 * @return the nonIdenticalAvgs
	 */
	public final double[] getNonIdenticalAvgs() {
		return nonIdenticalAvgs;
	}


	/**
	 * @return the nonInfomatives
	 */
	public final double[] getNonInfomatives() {
		return nonInfomatives;
	}


	/**
	 * @return the nonInfomativeAvgs
	 */
	public final double[] getNonInfomativeAvgs() {
		return nonInfomativeAvgs;
	}


	/**
	 * @return the nonMendelians
	 */
	public final double[] getNonMendelians() {
		return nonMendelians;
	}


	/**
	 * @return the nonMendelianAvgs
	 */
	public final double[] getNonMendelianAvgs() {
		return nonMendelianAvgs;
	}


	/**
	 * @return the dominantInheritanceStates
	 */
	public final QuartetInheritanceState[] getDominantInheritanceStates() {
		return dominantInheritanceStates;
	}
}
