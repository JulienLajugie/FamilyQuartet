package dataStructures;



/**
 * Continuous block on a chromosome having the same {@link QuartetInheritanceState} 
 * @author Julien Lajugie
 */
public class QuartetInheritanceStateBlock implements InheritanceStateBlock<QuartetInheritanceState> {
	
	private final QuartetInheritanceState 	blockState;				// inheritance state of the block
	private final String 					chromosome;				// chromosome of the block
	private final int 						startPosition;			// start position of the block
	private final int 						stopPosition;			// stop position of the block
	private int 							variantCount = 0;		// count of the variants in the block
	private int 							MIECount = 0;			// count of the variants with a mendelian inheritance error
	private int 							SCECount = 0;			// count of the variants with a state inconsistency error
	private int 							NICount = 0;			// count of the non informative variants


	/**
	 * Creates an instance of {@link QuartetInheritanceStateBlock}
	 * @param blockState {@link QuartetInheritanceState} of the block. Should be null if the block is partial
	 * @param chromosome chromosome of the block
	 * @param startPosition start position of the block
	 * @param stopPosition stop position of the block
	 */
	public QuartetInheritanceStateBlock(QuartetInheritanceState blockState, String chromosome, int startPosition, int stopPosition) {
		super();
		this.blockState = blockState;
		this.chromosome = chromosome;
		this.startPosition = startPosition;
		this.stopPosition = stopPosition;
	}


	/**
	 * Analyzes a variant to generate stats on the number of MIE / SCE / Not informative variant per block
	 * @param variant variant to analyze
	 */
	public void analyzeVariant(Variant variant) {
		// we add the variant to the list of variant of the block
		variantCount++;
		// we check if the variant is a MIE
		// if not check if it's a NI variant 
		// if not check if it's a SCE variant
		if (variant.isMIE()) {
			MIECount++;
		} else if (variant.isNotInformative()) {
			NICount++;
		} else if (isSCE(variant)) {
			SCECount++;
		}
	}


	/**
	 * @param variant a {@link Variant}
	 * @return true if the specified variant is a SCE (State Consistency Error)
	 */
	public boolean isSCE(Variant variant) {
		boolean isMIE = false;
		boolean isNI = false;
		boolean isSCE = true;
		for (QuartetInheritanceState currentInheritanceState: variant.getInheritanceStates()) {
			if (currentInheritanceState == QuartetInheritanceState.MIE) {
				isMIE = true;
			}
			if (currentInheritanceState == QuartetInheritanceState.NOT_INFORMATIVE) {
				isNI = true;
			}
			if (currentInheritanceState == blockState) {
				isSCE = false;
			}
		}
		return ((!isMIE) && (!isNI) && (isSCE));
	}


	/**
	 * @return the blockState
	 */
	public final QuartetInheritanceState getBlockState() {
		return blockState;
	}


	/**
	 * @return the chromosome
	 */
	public final String getChromosome() {
		return chromosome;
	}


	/**
	 * @return the startPosition
	 */
	public final int getStartPosition() {
		return startPosition;
	}


	/**
	 * @return the stopPosition
	 */
	public final int getStopPosition() {
		return stopPosition;
	}


	/**
	 * @return the number of variants in the block
	 */
	public final int getVariantCount() {
		return variantCount;
	}


	/**
	 * @return the MIECount
	 */
	public final int getMIECount() {
		return MIECount;
	}


	/**
	 * @return the SCECout
	 */
	public final int getSCECount() {
		return SCECount;
	}


	/**
	 * @return the NICount
	 */
	public int getNICount() {
		return NICount;
	}


	/**
	 * @return the percentage of MIE in the block
	 */
	public double getMIEPercentage() {
		if (variantCount == 0) {
			return 0;
		}
		return MIECount / (double) variantCount * 100d;
	}


	/**
	 * @return the percentage of SCE in the block
	 */
	public double getSCEPercentage() {
		if (variantCount == 0) {
			return 0;
		}
		return SCECount / (double) variantCount * 100d;
	}


	/**
	 * @return the percentage of non informative variants in the block
	 */
	public double getNIPercentage() {
		if (variantCount == 0) {
			return 0;
		}
		return NICount / (double) variantCount * 100d;
	}


	/**
	 * @return a string containing the statistics of the bloc
	 */
	public String getStatistics() {
		String statistics = chromosome + "\t";
		statistics += startPosition +"\t";
		statistics += stopPosition +"\t";
		statistics += blockState +"\t";
		statistics += variantCount +"\t";
		statistics += MIECount +"\t";
		statistics += getMIEPercentage() +"\t";
		statistics += SCECount +"\t";
		statistics += getSCEPercentage() +"\t";
		statistics += NICount +"\t";
		statistics += getNIPercentage();
		return statistics;
	}
}
