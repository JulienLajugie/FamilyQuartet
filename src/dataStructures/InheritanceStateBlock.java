package dataStructures;

import java.util.ArrayList;
import java.util.List;



/**
 * Continuous block on a chromosome having the same {@link InheritanceState} 
 * @author Julien Lajugie
 */
public class InheritanceStateBlock {

	private static final String STAT_HEADER = "chromosome\tstart\tstop\tstate\tvariant#\tMIE#\tMIE%\tSCE#\tSCE%\tNI#\tNI%"; // header corresponding the the statistics information line
	private final InheritanceState 	blockState;			// inheritance state of the block
	private final boolean			isPartial;			// true if the block is partial
	private final String 			chromosome;			// chromosome of the block
	private final int 				startPosition;		// start position of the block
	private final int 				stopPosition;		// stop position of the block
	private final List<Variant> 	variantList;		// list of the variants in the block
	private int 					variantCount = 0;	// count of the variants in the block
	private int 					MIECount = 0;		// count of the variants with a mendelian inheritance error
	private int 					SCECount = 0;		// count of the variants with a state inconsistency error
	private int 					NICount = 0;		// count of the non informative variants


	/**
	 * Creates an instance of {@link InheritanceStateBlock}
	 * @param blockState {@link InheritanceState} of the block. Should be null if the block is partial
	 * @param chromosome chromosome of the block
	 * @param startPosition start position of the block
	 * @param stopPosition stop position of the block
	 */
	public InheritanceStateBlock(InheritanceState blockState, String chromosome, int startPosition, int stopPosition) {
		super();
		this.blockState = blockState;
		this.isPartial = (blockState == null);
		this.chromosome = chromosome;
		this.startPosition = startPosition;
		this.stopPosition = stopPosition;
		this.variantList = new ArrayList<Variant>();
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
		boolean isMIE = false;
		boolean isNI = false;
		boolean isSCE = true;
		for (InheritanceState currentInheritanceState: variant.getInheritanceStates()) {
			if (currentInheritanceState == InheritanceState.MIE) {
				isMIE = true;
			}
			if (currentInheritanceState == InheritanceState.NOT_INFORMATIVE) {
				isNI = true;
			}
			if (isPartial || (currentInheritanceState == blockState)) {
				isSCE = false;
			}
		}
		if (isMIE) {
			MIECount++;
		} else if (isNI) {
			NICount++;
		} else if (isSCE) {
			SCECount++;
		}
	}	


	/**
	 * This methods prints in bgr format in the standard output the coordinates of the specified variant only if it's a SCE 
	 * @param variant {@link Variant}
	 */
	public void printSCEVariantBrgFormat(Variant variant) {
		if (isSCE(variant)) {
			System.out.println(variant.getChromosome() + "\t" + variant.getPosition() + "\t" + variant.getPosition() + "\t1");
		}
	}

	
	
	/**
	 * This methods prints in vcf format in the standard output the vcf info of the specified variant only if it's a SCE 
	 * @param variant a {@link Variant}
	 * @param vcfLine line from the vcf from which the variant was extracted
	 */
	public void printSCEVariantVcfFormat(Variant variant, String vcfLine) {
		if (isSCE(variant)) {
			System.out.println(vcfLine);
		}
	}
	
	
	/**
	 * @param variant a {@link Variant}
	 * @return true if the specified variant is a SCE (State Consistency Error)
	 */
	public boolean isSCE(Variant variant) {
		if (isPartial) {
			return false;
		}
		boolean isMIE = false;
		boolean isNI = false;
		boolean isSCE = true;
		for (InheritanceState currentInheritanceState: variant.getInheritanceStates()) {
			if (currentInheritanceState == InheritanceState.MIE) {
				isMIE = true;
			}
			if (currentInheritanceState == InheritanceState.NOT_INFORMATIVE) {
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
	public final InheritanceState getBlockState() {
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
	 * @return the variantList
	 */
	public final List<Variant> getVariantList() {
		return variantList;
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
	 * @return the isPartial
	 */
	public boolean isPartial() {
		return isPartial;
	}


	/**
	 * @return the header corresponding to the statistic line
	 */
	public static String getStatisticsHeader() {
		return STAT_HEADER;
	}


	/**
	 * @return a string containing the statistics of the bloc
	 */
	public String getStatistics() {
		String statistics = chromosome + "\t";
		statistics += startPosition +"\t";
		statistics += stopPosition +"\t";
		if (isPartial) {
			statistics += "partial\t";
		} else {
			statistics += blockState +"\t";
		}
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
