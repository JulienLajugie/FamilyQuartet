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
	private final List<Integer> positions;	// list containing the positions associated to the genotypes
	private String chromosome;				// chromosome of the series
	private boolean isAnalyzed;				// true if the genotypes list has been analyzed
	private int compatibleGenotypes;		// count of compatible genotype vectors
	private int incompatibleGenotypes;		// count of incompatible genotype vectors
	private Integer firstPosition;
	private String firstGenotype;



	/**
	 * Creates an instance of {@link PhasedGenotypesSeries}
	 */
	public PhasedGenotypesSeries() {
		super();
		genotypes = new ArrayList<String>();
		positions = new ArrayList<Integer>();
		compatibleGenotypes = 0;
		incompatibleGenotypes = 0;
		isAnalyzed = false;
		chromosome = null;
		firstPosition = null;
		firstGenotype = null;
	}


	/**
	 * Adds two genetic genotypes to the series of phased genotypes if both genotypes are phased and if they are heterozygous 
	 * @param genotype1 first genetic genotype to add
	 * @param genotype2 second genetic genotype to add
	 * @param chromosome chromosome of the genotypes
	 * @param position position of the genotypes
	 * @return GENOTYPES_ADDED if the genotypes were added,
	 * NOTHING_TO_ADD if the genotypes are not heterozygous,
	 * SERIES_FINISHED if at least one of the genotypes is not phased 
	 */
	public int add2GeneticGenotypes(String genotype1, String genotype2, String chromosome, int position) {
		// case where we start a new chromosome
		if (this.chromosome == null) {
			this.chromosome = chromosome;
		} else if (!this.chromosome.equals(chromosome)) {
			return SERIES_FINISHED;
		}
		// we just want to study heterozygous genotypes
		if (genotype1.charAt(0) == genotype1.charAt(2)) {
			return NOTHING_TO_ADD;
		}
		// case where at least one of the genotype is not phased
		if ((genotype1.charAt(1) != '|') || (genotype2.charAt(1) != '|')) {
			return SERIES_FINISHED;
		}
		positions.add(position);
		genotypes.add(genotype1 + genotype2);
		return GENOTYPES_ADDED;
	}


	/**
	 * Adds two genetic genotypes to the series of phased genotypes if both genotypes are phased and if they are heterozygous 
	 * @param geneticGenotype genetic genotype to add
	 * @param physicalGenotype physical genotype to add
	 * @param chromosome chromosome of the genotypes
	 * @param position position of the genotypes
	 * @return GENOTYPES_ADDED if the genotypes were added,
	 * NOTHING_TO_ADD if the genotypes are not heterozygous,
	 * SERIES_FINISHED if at least one of the genotypes is not phased 
	 */
	public int addGeneticPhysicalGenotypes(String geneticGenotype, String physicalGenotype, String chromosome, int position) {
		// case where we start a new chromosome
		if (this.chromosome == null) {
			this.chromosome = chromosome;
		} else if (!this.chromosome.equals(chromosome)) {
			if (firstPosition != null) {
				positions.add(firstPosition);
				genotypes.add(firstGenotype);
				firstPosition = null;
				firstGenotype = null;
			}
			return SERIES_FINISHED;
		}
		// we just want to study heterozygous genotypes
		if (geneticGenotype.charAt(0) == geneticGenotype.charAt(2)) {
			return NOTHING_TO_ADD;
		}
		// case where the genetic genotype is not phased
		if (geneticGenotype.charAt(1) != '|') {
			if (firstPosition != null) {
				positions.add(firstPosition);
				genotypes.add(firstGenotype);
				firstPosition = null;
				firstGenotype = null;
			}
			return SERIES_FINISHED;
		} else if (physicalGenotype.charAt(1) != '|') {
			// case where the physical genotype is not phased
			firstPosition = position;
			firstGenotype = geneticGenotype + physicalGenotype;
			return SERIES_FINISHED;
		}
		// case where both genotypes are phased
		if (firstPosition != null) {
			positions.add(firstPosition);
			genotypes.add(firstGenotype);
			firstPosition = null;
			firstGenotype = null;
		}
		positions.add(position);
		genotypes.add(geneticGenotype + physicalGenotype);
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
		positions.clear();
		chromosome = null;
	}


	/**
	 * Computes the number of elements with a compatible phasing and the number of element with an incompatible phasing 
	 */
	private void analyzeGenotypes() {
		if (genotypes.size() > 1) {
			int compatibleGenotypesTmp = 0;
			int incompatibleGenotypesTmp = 0;
			int homozygousVectorsCount = 0;
			for (String currentGenotype: genotypes) {
				if (currentGenotype.charAt(0) == currentGenotype.charAt(2)) {
					homozygousVectorsCount++;
				} else {
					if((currentGenotype.charAt(0) == currentGenotype.charAt(3)) && (currentGenotype.charAt(2) == currentGenotype.charAt(5))) {
						compatibleGenotypesTmp++;
					} else {
						incompatibleGenotypesTmp++;
					}
				}
			}
			compatibleGenotypes = Math.max(compatibleGenotypesTmp, incompatibleGenotypesTmp) + homozygousVectorsCount;
			incompatibleGenotypes = Math.min(compatibleGenotypesTmp, incompatibleGenotypesTmp);
		}
		isAnalyzed = true;
	}


	/**
	 * Prints the result of the phasing in a bedgraph format.
	 * A line with a score of 1 correspond to a line with a compatible phasing
	 * A line with a score of -1 correspond to a line with an incompatible phasing
	 */
	public void printResultPhasingBgr() {
		if (!isAnalyzed) {
			analyzeGenotypes();
		}
		List<Integer> compatibleList = new ArrayList<Integer>();
		List<Integer> incompatibleList = new ArrayList<Integer>();
		//if (genotypes.size() > 1) {
			for (int i = 0; i < genotypes.size(); i++) {
				String currentGenotype = genotypes.get(i);
				Integer currentPosition = positions.get(i);

				if((currentGenotype.charAt(0) == currentGenotype.charAt(3)) && (currentGenotype.charAt(2) == currentGenotype.charAt(5))) {
					compatibleList.add(currentPosition);
				} else {
					incompatibleList.add(currentPosition);
				}
			}
			// if there are more incompatible than compatible vectors we invert the lists
			if (incompatibleList.size() > compatibleList.size()) {
				List<Integer> listTmp = compatibleList;
				compatibleList = incompatibleList;
				incompatibleList = listTmp;
			}
			// print on the bedgraph format with a score of 1 for compatible vector or -1 for incompatible
			for (int i = 0; i < positions.size(); i++) {				
				int currentPosition = positions.get(i);
				int score = 0;
				if (compatibleList.contains(currentPosition)) {
					score = 1;
				} else {
					score = -1;
				}
				// the first element has a double score
				if (i == 0) {
					score *= 2;
				}
				// the first element has a triple score
				if (i == positions.size() - 1) {
					score *= 3;
				}
				System.out.println(chromosome + '\t' + currentPosition + '\t' + (currentPosition + 1) + "\t" + score);
			}
		//}
	}


	/**
	 * @return a segmental duplication block with the start and the stop position of the series
	 */
	public SegmentalDuplication getBlock() {
		if (!positions.isEmpty()) {
			int blockStart = positions.get(0);
			int blockStop = positions.get(positions.size() - 1);
			return new SegmentalDuplication(blockStart, blockStop);
		}
		return null;
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
