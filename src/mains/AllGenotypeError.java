package mains;

import dataStructures.AlleleType;
import dataStructures.CrossTriosInheritanceState;
import dataStructures.QuartetMember;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;


/**
 * @author Julien Lajugie
 */
public class AllGenotypeError {


	public static void main(String[] args) {
		allGenotypeError();
	}


	public static void allGenotypeError() {
		CrossTriosInheritanceState identicalState = new CrossTriosInheritanceState(TrioInheritanceState.IDENTICAL, TrioInheritanceState.IDENTICAL);
		CrossTriosInheritanceState nonIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.NON_IDENTICAL, TrioInheritanceState.NON_IDENTICAL);
		CrossTriosInheritanceState paternalIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.IDENTICAL, TrioInheritanceState.NON_IDENTICAL);
		CrossTriosInheritanceState maternalIdenticalState = new CrossTriosInheritanceState(TrioInheritanceState.NON_IDENTICAL, TrioInheritanceState.IDENTICAL);
		CrossTriosInheritanceState[] allStates = {identicalState, nonIdenticalState, paternalIdenticalState, maternalIdenticalState};

		int detectedAsMIECount = 0;
		int detectedAsSCECount = 0;
		int notDetectedCount = 0;

		for (int currentGenotype = 0; currentGenotype < 256; currentGenotype++) {
			Variant currentVariant = new Variant(currentGenotype);
			//if (!isAltRefGenotype(currentVariant)) {
				if (!currentVariant.isMIE()) {
					for (CrossTriosInheritanceState currentState: allStates) {
						if (!currentVariant.isSCE(currentState)) {
							int[] errorFromCurrentGenotype = generateErrorGenotypes(currentGenotype);
							for (int currentErrorGenotype: errorFromCurrentGenotype) {
								Variant currentErrorVariant = new Variant(currentErrorGenotype);
								if (currentErrorVariant.isMIE()) {
									detectedAsMIECount++;
								} else if (currentErrorVariant.isSCE(currentState)) {
									detectedAsSCECount++;
								} else {
									notDetectedCount++;
								}
							}
						}
					}
				}
			//}
		}
		int errorCount = notDetectedCount + detectedAsMIECount + detectedAsSCECount;
		double detectedAsMIEPercentage = detectedAsMIECount / (double) errorCount *100d;
		double detectedAsSCEPercentage = detectedAsSCECount / (double) errorCount *100d;
		double notDetectedPercentage = notDetectedCount / (double) errorCount *100d;

		System.out.println("Error#:\t" + errorCount);
		System.out.println("Error detected as MIE#:\t" + detectedAsMIECount + "\t%:\t" + detectedAsMIEPercentage);
		System.out.println("Error detected as SCE#:\t" + detectedAsSCECount + "\t%:\t" + detectedAsSCEPercentage);
		System.out.println("Error Not detected#:\t" + notDetectedCount + "\t%:\t" + notDetectedPercentage);
	}


	private static int[] generateErrorGenotypes(int currentGenotype) {
		int[] errorGenotypes = new int[8];		
		for (int i = 0; i < 8; i++) {
			errorGenotypes[i] = (currentGenotype ^ (0x1 << i));
		}
		return errorGenotypes;
	}

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
