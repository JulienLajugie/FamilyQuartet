package dataStructures;


/**
 * Continuous block on a chromosome having the same {@link QuartetInheritanceState}
 * @author Julien Lajugie
 * @param <T> type of inheritance state block (must implement the interface {@link InheritanceState})
 */
public interface InheritanceStateBlock<T extends InheritanceState> {
	
	/**
	 * Header corresponding the the statistics information line
	 */
	public static final String STAT_HEADER = "chromosome\tstart\tstop\tstate\tvariant#\tMIE#\tMIE%\tSCE#\tSCE%\tNI#\tNI%";

	
	/**
	 * Analyzes a variant to generate stats on the number of MIE / SCE / Not informative variant per block
	 * @param variant variant to analyze
	 */
	public void analyzeVariant(Variant variant);


	/**
	 * @param variant a {@link Variant}
	 * @return true if the specified variant is a SCE (State Consistency Error)
	 */
	public boolean isSCE(Variant variant);


	/**
	 * @return the blockState
	 */
	public T getBlockState();


	/**
	 * @return the chromosome
	 */
	public String getChromosome();


	/**
	 * @return the startPosition
	 */
	public int getStartPosition();


	/**
	 * @return the stopPosition
	 */
	public int getStopPosition();


	/**
	 * @return the number of variants in the block
	 */
	public int getVariantCount();


	/**
	 * @return the MIECount
	 */
	public int getMIECount();


	/**
	 * @return the SCECout
	 */
	public int getSCECount();


	/**
	 * @return the NICount
	 */
	public int getNICount();


	/**
	 * @return the percentage of MIE in the block
	 */
	public double getMIEPercentage();


	/**
	 * @return the percentage of SCE in the block
	 */
	public double getSCEPercentage();


	/**
	 * @return the percentage of non informative variants in the block
	 */
	public double getNIPercentage();


	/**
	 * @return a string containing the statistics of the bloc
	 */
	public String getStatistics();
}
