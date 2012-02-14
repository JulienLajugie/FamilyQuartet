package dataStructures;

import java.util.ArrayList;
import java.util.List;


/**
 * This class offers tools to compare the result of the phasing made by two different methods (eg: read backed phasing and haplotyping phasing).
 * @author Julien Lajugie
 */
public class PhasedGenotypesSeries {

	/**
	 * The genotypes were added
	 */
	public static final int GENOTYPES_ADDED = 0;

	/**
	 * At least one of the genotypes is not phased (so the phased genotype series ends)
	 */
	public static final int SERIES_FINISHED = 1;

	/**
	 * The genotypes are not heterozygous (no phasing to analyze)
	 */
	public static final int NOTHING_TO_ADD = 2;

	private final List<String> genotypes; 	// list containing the genotype vectors (each element is the concatenation of the 2 genotypes to compare)
	private String lastGenotypes;			// value of the last genotype vectors not added (because they were not phased) 
	private boolean isAnalyzed;				// true if the genotypes list has been analyzed
	private int compatibleGenotypes;		// count of compatible genotype vectors
	private int incompatibleGenotypes;		// count of incompatible genotype vectors


	/**
	 * Creates an instance of {@link PhasedGenotypesSeries}
	 */
	public PhasedGenotypesSeries() {
		super();
		genotypes = new ArrayList<String>();
		compatibleGenotypes = 0;
		incompatibleGenotypes = 0;
		isAnalyzed = false;
		lastGenotypes = null;
	}


	/**
	 * Adds two genotypes to the series of phased genotypes if both genotypes are phased and if they are heterozygous 
	 * @param genotype1 first genotype to add
	 * @param genotype2 second genotype to add
	 * @return GENOTYPES_ADDED if the genotypes were added,
	 * NOTHING_TO_ADD if the genotypes are not heterozygous,
	 * SERIES_FINISHED if at least one of the genotypes is not phased 
	 */
	public int addGenotypes(String genotype1, String genotype2) {
		// we just want to study heterozygous genotypes
		if (genotype1.charAt(0) == genotype1.charAt(2)) {
			lastGenotypes = null;
			return NOTHING_TO_ADD;
		}
		// case where at least one of the genotype is not phased
		if ((genotype1.charAt(1) != '|') || (genotype2.charAt(1) != '|')) {
			lastGenotypes = genotype1 + genotype2;
			return SERIES_FINISHED;
		}
		// case where both genotypes are phased
		if ((genotypes.size() == 0) && (lastGenotypes != null)) {
			genotypes.add(lastGenotypes);	
		}
		genotypes.add(genotype1 + genotype2);
		return GENOTYPES_ADDED;
	}


	/**
	 * Resets the genotypes series.  Removes the stored genotypes and set the statistics values to 0
	 */
	public void reset() {
		compatibleGenotypes = 0;
		incompatibleGenotypes = 0;
		isAnalyzed = false;
		genotypes.clear();
	}


	/**
	 * Computes the number of elements with a compatible phasing and the number of element with an incompatible phasing 
	 */
	private void analyzeGenotypes() {
		if (genotypes.size() > 1) {
			int compatibleGenotypesTmp = 0;
			int incompatibleGenotypesTmp = 0;
			for (String currentGenotype: genotypes) {
				if((currentGenotype.charAt(0) == currentGenotype.charAt(3)) && (currentGenotype.charAt(2) == currentGenotype.charAt(5))) {
					compatibleGenotypesTmp++;
				} else {
					incompatibleGenotypesTmp++;
				}
			}
			compatibleGenotypes = Math.max(compatibleGenotypesTmp, incompatibleGenotypesTmp);
			incompatibleGenotypes = Math.min(compatibleGenotypesTmp, incompatibleGenotypesTmp);
		}
		isAnalyzed = true;
	}


	/**
	 * @return the count of compatible genotypes
	 */
	public int getCompatibleGenotypes() {
		if (!isAnalyzed) {
			analyzeGenotypes();
		}
		return compatibleGenotypes;
	}


	/**
	 * @return the count of incompatible genotypes
	 */
	public int getIncompatibleGenotypes() {
		if (!isAnalyzed) {
			analyzeGenotypes();
		}
		return incompatibleGenotypes;
	}


	/**
	 * Prints the series
	 */
	public void print() {
		for (String currentGenotype: genotypes) {
			System.out.println(currentGenotype);
		}
		System.out.println("Compatible vector count:" + getCompatibleGenotypes() + ", incompatible vector count:" + getIncompatibleGenotypes());
	}
}
