package mains;

import dataStructures.AlleleType;
import dataStructures.CrossTriosInheritanceState;
import dataStructures.QuartetMember;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;


/**
 * Generates all the possible genotypes
 * For each genotype generates all the possible single errors
 * For each inheritance state where the genotype is not a SCE nor a MIE, 
 * checks if the error is detected as a SCE or a MIE or if the error is not detected
 * @author Julien Lajugie
 */
public class AllGenotypeError {

	/**
	 * Main method
	 * @param args
	 */
	public static void main(String[] args) {
		allGenotypeError();
	}


	/**
	 *Generates all the possible genotypes
	 * For each genotype generates all the possible single errors
	 * For each inheritance state where the genotype is not a SCE nor a MIE, 
	 * checks if the error is detected as a SCE or a MIE or if the error is not detected 
	 */
	public static void allGenotypeError() {
		CrossTriosInheritanceState identicalState = new CrossTriosInheritanceState(TrioInheritanceState.IDENTICAL, TrioInheritanceState.IDENTICAL);
		CrossTriosInheritanceState nonIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.NON_IDENTICAL, TrioInheritanceState.NON_IDENTICAL);
		CrossTriosInheritanceState paternalIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.IDENTICAL, TrioInheritanceState.NON_IDENTICAL);
		CrossTriosInheritanceState maternalIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.NON_IDENTICAL, TrioInheritanceState.IDENTICAL);
		CrossTriosInheritanceState[] allStates = {identicalState, nonIdenticalState, paternalIdenticalState, maternalIdenticalState};

		int detectedAsMIECount = 0;
		int detectedAsSCECount = 0;
		int notDetectedCount = 0;

		int parentError = 0;
		int kidsError = 0;
		int test = 1;
		for (int currentGenotype = 0; currentGenotype < 256; currentGenotype++) {
			Variant currentVariant = new Variant(currentGenotype);
			// we want to eliminate the 1/0 genotypes since they are equivalent to the 0/1 ones
			if (!isAltRefGenotype(currentVariant)) {

				if (!currentVariant.isMIE()) {
					for (CrossTriosInheritanceState currentState: allStates) {
						
					//CrossTriosInheritanceState currentState = paternalIdenticalState;
						if (!currentVariant.isSCE(currentState)) {
							System.out.println(test++);

							int[] errorFromCurrentGenotype = generateErrorGenotypes(currentGenotype);
							int indexError = 0;
							for (int currentErrorGenotype: errorFromCurrentGenotype) {
								Variant currentErrorVariant = new Variant(currentErrorGenotype);
		
								if (currentErrorVariant.isMIE()) {

									detectedAsMIECount++;
								} else if (currentErrorVariant.isSCE(currentState)) {

									detectedAsSCECount++;
								} else {
									if (indexError >= 4) {
										parentError++;
									} else {
										kidsError++;
									}
									notDetectedCount++;
								}
								indexError++;
							}
						}
					}
				}
			}
		}
		int errorCount = notDetectedCount + detectedAsMIECount + detectedAsSCECount;
		double detectedAsMIEPercentage = detectedAsMIECount / (double) errorCount *100d;
		double detectedAsSCEPercentage = detectedAsSCECount / (double) errorCount *100d;
		double notDetectedPercentage = notDetectedCount / (double) errorCount *100d;

		System.out.println("Error#:\t" + errorCount);
		System.out.println("Error detected as MIE#:\t" + detectedAsMIECount + "\t%:\t" + detectedAsMIEPercentage);
		System.out.println("Error detected as SCE#:\t" + detectedAsSCECount + "\t%:\t" + detectedAsSCEPercentage);
		System.out.println("Error Not detected#:\t" + notDetectedCount + "\t%:\t" + notDetectedPercentage);
		
		System.out.println("parent Errors = " + parentError);
		System.out.println("kids Errors = " + kidsError);
	}


	/**
	 * @param currentGenotype initial genotype
	 * @return an array with all the possible error;
	 * 8 elements, an error can occure at any of the 8 allele
	 */
	private static int[] generateErrorGenotypes(int currentGenotype) {
		int[] errorGenotypes = new int[8];
		for (int i = 0; i < 8; i++) {
			errorGenotypes[i] = (currentGenotype ^ (0x1 << i));
		}
		return errorGenotypes;
	}


	/**
	 * @param variant
	 * @return true if one of the familly member has a 1/0 genotype
	 */
	private static boolean isAltRefGenotype(Variant variant) {
		for (QuartetMember currentMember: QuartetMember.values()) {
			AlleleType[] alleles = variant.getAlleles(currentMember);
			if ((alleles[0] == AlleleType.ALTERNATIVE_ALLELE) && (alleles[1] == AlleleType.REFERENCE_ALLELE)) {
				return true;
			}
		}
		return false;
	}
}
