package utils;

import java.util.HashMap;
import java.util.Map;

import dataStructures.InheritanceState;



/**
 * This class allows to retrieve the inheritance states from the genotype pattern as describe in 
 * the figure S2A of the supplementary data of the following article: 
 * "Analysis of Genetic Inheritance in a Family Quartet by Whole-Genome Sequencing", Roach et Al
 * @author Julien Lajugie
 */
public class PatternToInheritanceStates {
	
	/**
	 *  hash map with the value for the conversion genotype pattern to inheritance states
	 */
	private static final Map<String, InheritanceState[]> MAP = new HashMap<String, InheritanceState[]>(); 

	// populate the map
	static {
		/*	Ten possible family genotype patterns for biallelic
		positions are consistent with one or two inheritance states. For each pattern, 
		the parental genotypes are shown on top (father to the left, mother to the right).
		The most frequent allele in the parents is denoted by “a”; in case of equal frequency, “a”
		denotes the most frequent allele in the children. Thus, “aa+ab” means the father is
		homozygous for the most frequent allele, while the mother is heterozygous. A “/” symbol
		is used to indicate that the order is not important.*/
		MAP.put("ab/ab;aa/aa", new InheritanceState[]{InheritanceState.IDENTICAL});
		MAP.put("ab+aa;aa/aa", new InheritanceState[]{InheritanceState.IDENTICAL, InheritanceState.PATERNAL});
		MAP.put("aa+ab;aa/aa", new InheritanceState[]{InheritanceState.IDENTICAL, InheritanceState.MATERNAL});
		MAP.put("ab+aa;ab/ab", new InheritanceState[]{InheritanceState.IDENTICAL, InheritanceState.PATERNAL});
		MAP.put("aa+ab;ab/ab", new InheritanceState[]{InheritanceState.IDENTICAL, InheritanceState.MATERNAL});
		MAP.put("ab/ab;aa/ab", new InheritanceState[]{InheritanceState.PATERNAL, InheritanceState.MATERNAL});
		MAP.put("ab/ab;ab/ab", new InheritanceState[]{InheritanceState.IDENTICAL, InheritanceState.NON_IDENTICAL});
		MAP.put("aa+ab;aa/ab", new InheritanceState[]{InheritanceState.NON_IDENTICAL, InheritanceState.PATERNAL});
		MAP.put("ab+aa;aa/ab", new InheritanceState[]{InheritanceState.NON_IDENTICAL, InheritanceState.MATERNAL});
		MAP.put("ab/ab;aa/bb", new InheritanceState[]{InheritanceState.NON_IDENTICAL});
		
		
		/* One biallelic family genotype pattern, and the single monoallelic pattern, 
		 * are consistent with all inheritance states, and therefore uninformative. */
		MAP.put("aa/bb;ab/ab", new InheritanceState[]{InheritanceState.NOT_INFORMATIVE});
		MAP.put("aa/aa;aa/aa", new InheritanceState[]{InheritanceState.NOT_INFORMATIVE});
		
		/* Five genotype patterns are Mendelian Inheritance Errors (MIE), with a novel allele in a child. */
		MAP.put("aa/aa;aa/ab", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/aa;aa/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/aa;ab/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/aa;bb/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/aa;ab/ab", new InheritanceState[]{InheritanceState.MIE});
		
		/* Six genotype patterns would require that both alleles observed in a child derive from the same parent. */
		MAP.put("aa/bb;aa/aa", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/bb;aa/ab", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/bb;aa/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/ab;aa/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/ab;ab/bb", new InheritanceState[]{InheritanceState.MIE});
		MAP.put("aa/ab;bb/bb", new InheritanceState[]{InheritanceState.MIE});
	}
	
	
	/**
	 * @param genotypePattern a genotype pattern (ie aa+ab;aa/aa)
	 * @return the genotype state associated to the specified pattern
	 */
	public static InheritanceState[] getInheritanceStates(String genotypePattern) {
		return MAP.get(genotypePattern);
	}
}
