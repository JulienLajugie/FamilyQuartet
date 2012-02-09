package dataStructures;

import dataStructures.InheritanceState;
import exceptions.InvalidVCFLineException;
import utils.PatternToInheritanceStates;

/**
 * This class represents a variant from a VCF file
 * @author Julien Lajugie
 */
public class Variant {

	/**
	 * Filter using the GATK hard filtering recommendation from the best practice V3 if true
	 */
	private static final boolean USE_GATK_HARD_FILTERING = false;
	/**
	 * Filter using the some of the individual quality scores if true
	 */
	private static final boolean USE_INDIVIDUALS_SCORE_FILTERING = false;
	/**
	 * Filter using the min of the individual PL fields
	 */
	private static final boolean USE_INDIVIDUALS_PL_FILTERING = true;
	/**
	 * Filter removing all the variant with the filter field different from "PASS"
	 */
	private static final boolean USE_FILTER_FIELD_FILTERING = true;
	/**
	 * Filter excluding the fully heterozygous genotypes
	 */
	private static final boolean USE_FILTER_HETEROZYGOUS = false;
	/**
	 * Filter excluding the 3/4 heterozygous genotypes
	 */
	private static final boolean USE_FILTER_3QUATER_HETEROZYGOUS = false;
	
	
	private final String 				chromosome;				// chromosome of the variant		
	private final int 					position;				// position of the variant
	private final String 				referenceAllele;		// reference allele of the variant
	private final String 				alternatievAllele;		// alternative allele of the variant
	private final AlleleType[] 			fatherAlleles;			// alleles of the father
	private final AlleleType[] 			motherAlleles;			// alleles of the mother
	private final AlleleType[] 			kid1Alleles;			// alleles of the 1st kid
	private final AlleleType[] 			kid2Alleles;			// alleles of the 2nd kid
	private final String				genotypePattern;		// genotype pattern of the variant for the familly quartet
	private final InheritanceState[]	inheritanceStates;		// inheritance state of the variant


	/**
	 * Creates an instance of {@link Variant}
	 * @param VCFLine
	 * @throws InvalidVCFLineException when the VCF line is not valid
	 */
	@SuppressWarnings("unused")
	public Variant(String VCFLine) throws InvalidVCFLineException {
		String[] splitLine = VCFLine.split("\t"); // VCF fields are tab delimited
		// filter using the filter field
		if (USE_FILTER_FIELD_FILTERING) {
			filterFieldFiltering(splitLine[6].trim());
		}
		// filter using the GATK hard filtering recommendation from the best practice V3
		if (USE_GATK_HARD_FILTERING) {
			GATKFilter(splitLine[7].trim());
		}
		// we retrive the info about the chromosome, the position, the reference allele and the alternative allele
		chromosome = splitLine[0].trim();
		position = Integer.parseInt(splitLine[1].trim());
		referenceAllele = splitLine[3].trim();
		alternatievAllele = splitLine[4].trim();
		// filter using the some of the individual quality scores
		if (USE_INDIVIDUALS_SCORE_FILTERING) {
			double score = stringToQualityScore(splitLine[9].trim());
			score += stringToQualityScore(splitLine[10].trim());
			score += stringToQualityScore(splitLine[11].trim());
			score += stringToQualityScore(splitLine[12].trim());
			if (score < 200d) {
				throw new InvalidVCFLineException();
			}
			/*double indScore = Math.min(stringToQualityScore(splitLine[9].trim()), stringToQualityScore(splitLine[10].trim()));
			indScore = Math.min(indScore, stringToQualityScore(splitLine[11].trim()));
			indScore = Math.min(indScore, stringToQualityScore(splitLine[12].trim()));
			if (indScore < 25d) {
				throw new InvalidVCFLineException();
			}*/
		}
		// Filter using the min of the individual PL fields
		if (USE_INDIVIDUALS_PL_FILTERING) {
			int plScore = Math.min(formatFieldToPL(splitLine[9].trim()), formatFieldToPL(splitLine[10].trim()));
			plScore = Math.min(plScore, formatFieldToPL(splitLine[11].trim()));
			plScore = Math.min(plScore, formatFieldToPL(splitLine[12].trim()));
			if (plScore < 20) {
				throw new InvalidVCFLineException();
			}
		}
		// extract the allele informations
		fatherAlleles = stringToAlleleTypes(splitLine[9].trim());
		motherAlleles = stringToAlleleTypes(splitLine[10].trim());
		kid1Alleles = stringToAlleleTypes(splitLine[11].trim());
		kid2Alleles = stringToAlleleTypes(splitLine[12].trim());		
		// compute the genotype pattern
		genotypePattern = computeGenotypePattern();		
		// compute the inheritance states
		inheritanceStates = PatternToInheritanceStates.getInheritanceStates(genotypePattern);
		// exclude the fully heterozygote vectors if the filter is set to true
		if (USE_FILTER_HETEROZYGOUS && genotypePattern.equals("ab/ab;ab/ab")) {
			throw new InvalidVCFLineException(); 
		}
		// exclude the 3/4 heterozygote vectors if the filter is set to true
		if (USE_FILTER_3QUATER_HETEROZYGOUS && 
				((genotypePattern.equals("ab+aa;ab/ab") || genotypePattern.equals("aa+ab;ab/ab") || genotypePattern.equals("ab/ab;aa/ab")))) {
			throw new InvalidVCFLineException();
		}
	}


	/**
	 * 
	 * @param formatField format field following the GT:AD:DP:GQ:PL format
	 * @return the minimum PL score (ie the max probability that the genotype is not the one returned by the genotyper)
	 * @throws InvalidVCFLineException if there is more than 1 alternative
	 */
	private int formatFieldToPL(String formatField) throws InvalidVCFLineException {
		String[] splitFormatField = formatField.split(":");
		// the following happens when we have a ./. variant
		if (splitFormatField.length != 5) {
			throw new InvalidVCFLineException();
		}
		String genotype = splitFormatField[0].trim();
		String[] plScores = splitFormatField[4].trim().split(",");
		int refRefScore = Integer.parseInt(plScores[0].trim());
		int refAltScore = Integer.parseInt(plScores[1].trim());
		int altAltScore = Integer.parseInt(plScores[2].trim());
		if (genotype.equals("0/0")) {
			return Math.min(refAltScore, altAltScore);
		}
		if (genotype.equals("0/1")) {
			return Math.min(refRefScore, altAltScore);
		}
		if (genotype.equals("1/1")) {
			return Math.min(refRefScore, refAltScore);
		}
		throw new InvalidVCFLineException();
	}
	
	
	/**
	 * Throws an exception if the specified string is different from "PASS"
	 * @param filterField
	 * @throws InvalidVCFLineException
	 */
	private void filterFieldFiltering(String filterField) throws InvalidVCFLineException {
		if ((!filterField.equalsIgnoreCase("PASS"))) {
			// && (!filterField.equalsIgnoreCase("TruthSensitivityTranche99.00to99.90"))) {
			throw new InvalidVCFLineException(); 
		}
	}


	/**
	 * Hard filtering as recommended in the GATK best practice V3
	 * @param infoField info field of the VCF line
	 * @throws InvalidVCFLineException
	 */
	private void GATKFilter(String infoField) throws InvalidVCFLineException {
		// filter on the the QD field
		int QDIndex = infoField.indexOf("QD=");
		if (QDIndex == -1) {
			throw new InvalidVCFLineException(); 
		} else {
			String QDStr = infoField.substring(QDIndex + 3);
			int indexSemicolon = QDStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFLineException();
			} else {
				QDStr = QDStr.substring(0, indexSemicolon);
				double QD = Double.parseDouble(QDStr);
				//System.out.println("QD=" + QD);
				if (QD < 8.0) {
					throw new InvalidVCFLineException();
				}
			}			
		}
		// filter on the the HRun field
		int HRunIndex = infoField.indexOf("HRun=");
		if (HRunIndex == -1) {
			throw new InvalidVCFLineException(); 
		} else {
			String HRunStr = infoField.substring(HRunIndex + 5);
			int indexSemicolon = HRunStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFLineException();
			} else {
				HRunStr = HRunStr.substring(0, indexSemicolon);
				double HRun = Double.parseDouble(HRunStr);
				//System.out.println("HRun=" + HRun);
				if (HRun > 5) {
					throw new InvalidVCFLineException();
				}
			}			
		}
		// filter on the the FS field
		int FSIndex = infoField.indexOf("FS=");
		if (FSIndex == -1) {
			throw new InvalidVCFLineException(); 
		} else {
			String FSStr = infoField.substring(FSIndex + 3);
			int indexSemicolon = FSStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFLineException();
			} else {
				FSStr = FSStr.substring(0, indexSemicolon);
				double FS = Double.parseDouble(FSStr);
				//System.out.println("FS=" + FS);	
				if (FS > 200) {
					throw new InvalidVCFLineException();
				}			
			}
		}
	}


	/**
	 * @param genotypeInfo genotype information field of a VCF line
	 * @return the allele types of a sample
	 * @throws InvalidVCFLineException when the genotype doesn't contains information about the allele
	 */
	private AlleleType[] stringToAlleleTypes(String genotypeInfo) throws InvalidVCFLineException {
		AlleleType[] resultAlleles = new AlleleType[2];
		String[] splitGenotypeInfo = genotypeInfo.split(":"); // the genotype info field is colon-separated
		String genotype = splitGenotypeInfo[0];	// the genotype is in the first field
		// the genotype is coded like follow 0/0, 0/1 (where 0 is the reference allele and 1 is the alternative allele)
		if (genotype.charAt(0) == '0') {
			resultAlleles[0] = AlleleType.REFERENCE_ALLELE;
		} else if (genotype.charAt(0) == '1') {
			resultAlleles[0] = AlleleType.ALTERNATIVE_ALLELE;
		} else {
			throw new InvalidVCFLineException();
		}
		if(genotype.charAt(2) == '0') {
			resultAlleles[1] = AlleleType.REFERENCE_ALLELE;
		} else if (genotype.charAt(2) == '1') {
			resultAlleles[1] = AlleleType.ALTERNATIVE_ALLELE;
		} else {
			throw new InvalidVCFLineException();
		}
		//System.out.print(genotype + "\t");
		return resultAlleles;
	}


	/**
	 * @param genotypeInfo Genotype Info field of the VCF
	 * @return the quality score from the field Genotype Info
	 */
	private double stringToQualityScore(String genotypeInfo) {
		String[] splitGenotypeInfo = genotypeInfo.split(":");
		if (splitGenotypeInfo.length > 3) {
			String scoreStr = splitGenotypeInfo[3];
			try {
				return Double.parseDouble(scoreStr);
			} catch (NumberFormatException e) {
				return 0;
			}
		} else {
			return 0.0;
		}
	}


	/**
	 * @return the genotype pattern of the quartet
	 */
	private String computeGenotypePattern() {
		AlleleType mostFrequentAllele = getMostFrequentAllele();
		String fatherPattern = getSamplePattern(fatherAlleles, mostFrequentAllele);
		String motherPattern = getSamplePattern(motherAlleles, mostFrequentAllele);
		String kid1Pattern = getSamplePattern(kid1Alleles, mostFrequentAllele);
		String kid2Pattern = getSamplePattern(kid2Alleles, mostFrequentAllele);
		String kidPattern = getKidPattern(kid1Pattern, kid2Pattern);
		String parentPattern = getParentPattern(fatherPattern, motherPattern, kidPattern);
		return parentPattern + ";" + kidPattern;
	}


	/**
	 * The most frequent allele in the parents is denoted by “a”; in case of equal frequency, “a”
	 * denotes the most frequent allele in the children.
	 * @return the type of the most frequent allele (reference or alternate)
	 */
	private AlleleType getMostFrequentAllele() {
		int refCount = 0; // counter for the reference allele
		refCount = (fatherAlleles[0] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
		refCount = (fatherAlleles[1] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
		refCount = (motherAlleles[0] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
		refCount = (motherAlleles[1] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
		switch (refCount) {
		case 0:
		case 1:
			return AlleleType.ALTERNATIVE_ALLELE;
		case 2:
			refCount = (kid1Alleles[0] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
			refCount = (kid1Alleles[1] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
			refCount = (kid2Alleles[0] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
			refCount = (kid2Alleles[1] == AlleleType.REFERENCE_ALLELE) ? refCount + 1 : refCount;
			if (refCount < 4) {
				return AlleleType.ALTERNATIVE_ALLELE;
			} else {
				return AlleleType.REFERENCE_ALLELE;
			}
		case 3:
		case 4:
			return AlleleType.REFERENCE_ALLELE;
		}
		// should not happen
		return null;
	}


	/**
	 * @param alleleTypes the alleles type of an individual of the quartet 
	 * @param mostFrequentAllele the most frequent allele of the quartet as specified in the getMostFrequentAllele() function
	 * @return the genotype pattern for an individual of the quartet
	 */
	public String getSamplePattern (AlleleType[] alleleTypes, AlleleType mostFrequentAllele) {
		int mostFrequentCount = 0; // counter for the number of most frequent allele
		mostFrequentCount = (alleleTypes[0] == mostFrequentAllele) ? mostFrequentCount + 1 : mostFrequentCount;
		mostFrequentCount = (alleleTypes[1] == mostFrequentAllele) ? mostFrequentCount + 1 : mostFrequentCount;
		switch (mostFrequentCount) {
		case 0:
			return "bb";
		case 1:
			return "ab";
		case 2:
			return "aa";
		}
		// should not happen
		return null;
	}


	/**
	 * @param kid1Pattern genotype pattern of the first kid
	 * @param kid2Pattern genotype pattern of the second kid
	 * @return the children's pattern
	 */
	private String getKidPattern (String kid1Pattern, String kid2Pattern) {
		if (kid1Pattern.equals("aa")) {
			return kid1Pattern + "/" + kid2Pattern;
		}
		if (kid1Pattern.equals("bb")) {
			return kid2Pattern + "/" + kid1Pattern;
		}
		if (kid2Pattern.equals("aa")) {
			return kid2Pattern + "/" + kid1Pattern;
		}
		if (kid2Pattern.equals("bb")) {
			return kid1Pattern + "/" + kid2Pattern;
		}
		// last possible case
		return "ab/ab";		
	}


	/**
	 * @param fatherPattern pattern of the father
	 * @param motherPattern pattern of the mother
	 * @param kidPattern pattern of the children
	 * @return the parent's pattern
	 */
	private String getParentPattern (String fatherPattern, String motherPattern, String kidPattern) {
		if (fatherPattern.equals(motherPattern)) {
			return fatherPattern + "/" + motherPattern;
		}
		if (fatherPattern.equals("bb")) {
			return motherPattern + "/" + fatherPattern;
		}
		if (motherPattern.equals("bb")) {
			return fatherPattern + "/" + motherPattern;
		} 
		if (kidPattern.equals("aa/bb") || kidPattern.equals("ab/bb") || kidPattern.equals("bb/bb")) {
			return "aa/ab";
		}
		return fatherPattern + "+" + motherPattern;		
	}


	/**
	 * @return the chromosome of the variant
	 */
	public final String getChromosome() {
		return chromosome;
	}


	/**
	 * @return the position of the variant
	 */
	public final int getPosition() {
		return position;
	}


	/**
	 * @return the reference allele of the variant
	 */
	public final String getReferenceAllele() {
		return referenceAllele;
	}


	/**
	 * @return the alternatiev allele of the variant
	 */
	public final String getAlternatievAllele() {
		return alternatievAllele;
	}


	/**
	 * @return the father's alleles
	 */
	public final AlleleType[] getFatherAlleles() {
		return fatherAlleles;
	}


	/**
	 * @return the mother's alleles
	 */
	public final AlleleType[] getMotherAlleles() {
		return motherAlleles;
	}


	/**
	 * @return the 1st kid's alleles
	 */
	public final AlleleType[] getKid1Alleles() {
		return kid1Alleles;
	}


	/**
	 * @return the 2nd kid's alleles
	 */
	public final AlleleType[] getKid2Alleles() {
		return kid2Alleles;
	}


	/**
	 * @return the genotype pattern of the variant for the family quartet
	 */
	public final String getGenotypePattern() {
		return genotypePattern;
	}


	/**
	 * @return the inheritance states of the variant for the family quartet
	 */
	public final InheritanceState[] getInheritanceStates() {
		return inheritanceStates;
	}


	/**
	 * @return a string with the different field names of a variant
	 */
	public static String variantHeader() {
		return "chromosome\t" +
				"position\t" +
				"reference allele\t" +
				"alternative allele\t" +
				"father allele1\t" +
				"father allele2\t" +
				"mother allele1\t" +
				"mother allel2\t" +
				"kid1 allele1\t" +
				"kid1 allele2\t" +
				"kid2 allele1\t" +
				"kid2 allele2\t" +
				"genotype pattern\t" +
				"genotype state1\t" +
				"genotype state2";
	}

	
	@Override
	public String toString() {
		String variantString = "";
		variantString += chromosome;
		variantString += "\t";
		variantString += position;
		variantString += "\t";
		variantString += referenceAllele;
		variantString += "\t";
		variantString += alternatievAllele;
		variantString += "\t";
		variantString += fatherAlleles[0];
		variantString += "\t";
		variantString += fatherAlleles[1];
		variantString += "\t";
		variantString += motherAlleles[0];
		variantString += "\t";
		variantString += motherAlleles[1];
		variantString += "\t";
		variantString += kid1Alleles[0];
		variantString += "\t";
		variantString += kid1Alleles[1];
		variantString += "\t";
		variantString += kid2Alleles[0];
		variantString += "\t";
		variantString += kid2Alleles[1];
		variantString += "\t";
		variantString += genotypePattern;
		variantString += "\t";
		variantString += inheritanceStates[0];
		variantString += "\t";
		if (inheritanceStates.length > 1) {
			variantString += inheritanceStates[1];
		} else {
			variantString += "-";
		}
		return variantString;
	}
	
	
	/**
	 * This methods prints in bgr format in the standard output the coordinates of the specified variant
	 * if the variant states correspond to the specified states
	 * @param states {@link InheritanceState}
	 */
	public void printVariantBgrFormat(InheritanceState... states) {
		/*int statesFoundCount = 0;		
		for (InheritanceState currentVariantState: this.getInheritanceStates()) {
			for (InheritanceState currentInputState: states) {
				if (currentVariantState == currentInputState) {
					statesFoundCount++;
				}
			}
		}
		if (statesFoundCount == states.length) {*/
			System.out.println(this.getChromosome() + "\t" + this.getPosition() + "\t" + this.getPosition() + "\t1");
		//}
	}
}
