package dataStructures;

import dataStructures.QuartetInheritanceState;
import exceptions.FilteredVCFLineException;
import exceptions.InvalidVCFFieldException;
import exceptions.InvalidVCFLineException;
import exceptions.PartiallyCalledVariantException;
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
	 * Value of the PL filter. Set to null to disable
	 */
	public static final Integer INDIVIDUALS_PL_MIN_VALUE = null;
	/**
	 * Filter removing all the variant with the filter field different from "PASS"
	 */
	private static final VCFFilterField FILTER_FIELD_FILTERING = VCFFilterField.NONE;
	/**
	 * Filter excluding the fully heterozygous genotypes
	 */
	private static final boolean USE_FILTER_HETEROZYGOUS_FILTERING = false;
	/**
	 * Filter excluding the 3/4 heterozygous genotypes
	 */
	private static final boolean USE_FILTER_3QUATER_HETEROZYGOUS_FILTERING = false;
	/** 
	 * Only consider the genotypes with a score above this threshold as phased
	 */
	private static final double PHASING_QUALITY_FILTER_THRESHOLD = 0;
	/**
	 * Filters out MIE variants if set to true
	 */
	private static final boolean USE_MIE_FILTERING = false;

	private final String 					chromosome;					// chromosome of the variant		
	private final int 						position;					// position of the variant
	private final String 					referenceAllele;			// reference allele of the variant
	private final String 					alternatievAllele;			// alternative allele of the variant
	private final AlleleType[] 				fatherAlleles;				// alleles of the father
	private final AlleleType[] 				motherAlleles;				// alleles of the mother
	private final AlleleType[] 				kid1Alleles;				// alleles of the 1st kid
	private final AlleleType[] 				kid2Alleles;				// alleles of the 2nd kid
	private boolean							isFatherPhased;				// true if the father is phased
	private boolean							isMotherPhased;				// true if the mother is phased
	private boolean							isKid1Phased;				// true if kid 1
	private boolean							isKid2Phased;				// true if kid 2
	private final String					genotypePattern;			// genotype pattern of the variant for the familly quartet
	private final QuartetInheritanceState[]	quartetInheritanceStates;	// inheritance state of the variant
	private final int 						phasingQualityIndex;		// index of the phasing quality field


	/**
	 * Creates an instance of {@link Variant}
	 * @param VCFLine
	 * @throws InvalidVCFLineException when the VCF line is not valid
	 * @throws FilteredVCFLineException 
	 * @throws InvalidVCFFieldException 
	 * @throws PartiallyCalledVariantException 
	 */
	@SuppressWarnings("unused")
	public Variant(String VCFLine) throws InvalidVCFLineException, FilteredVCFLineException, InvalidVCFFieldException, PartiallyCalledVariantException {
		String[] splitLine = VCFLine.split("\t"); // VCF fields are tab delimited
		// filter using the filter field
		filterFieldFiltering(splitLine[6].trim());

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
				throw new FilteredVCFLineException("Individual score sum", Double.toString(score));
			}
		}
		// Filter using the min of the individual PL fields
		if (INDIVIDUALS_PL_MIN_VALUE != null) {
			int plScore = Math.min(genotypeFieldToPL(splitLine[9].trim()), genotypeFieldToPL(splitLine[10].trim()));
			plScore = Math.min(plScore, genotypeFieldToPL(splitLine[11].trim()));
			plScore = Math.min(plScore, genotypeFieldToPL(splitLine[12].trim()));
			if (plScore < INDIVIDUALS_PL_MIN_VALUE) {
				throw new FilteredVCFLineException("PL", Integer.toString(plScore));
			}
		}
		if (VCFLine.contains("PhasingInconsistent")) {
			phasingQualityIndex = -1;
			isFatherPhased = false;
			isMotherPhased = false;
			isKid1Phased = false;
			isKid2Phased = false;
		} else {
			phasingQualityIndex = getPhasingQualityFieldIndex(splitLine[8].trim());
			isFatherPhased = isGenotypePhased(splitLine[9].trim());
			isMotherPhased = isGenotypePhased(splitLine[10].trim());
			isKid1Phased = isGenotypePhased(splitLine[11].trim());
			isKid2Phased = isGenotypePhased(splitLine[12].trim());
		}
		// extract the allele informations
		fatherAlleles = stringToAlleleTypes(splitLine[9].trim());
		motherAlleles = stringToAlleleTypes(splitLine[10].trim());
		kid1Alleles = stringToAlleleTypes(splitLine[11].trim());
		kid2Alleles = stringToAlleleTypes(splitLine[12].trim());
		// compute the genotype pattern
		genotypePattern = computeGenotypePattern();		
		// compute the inheritance states
		quartetInheritanceStates = PatternToInheritanceStates.getInheritanceStates(genotypePattern);
		if (USE_MIE_FILTERING && (quartetInheritanceStates[0] == QuartetInheritanceState.MIE)) {
			throw new FilteredVCFLineException("MIE", "MIE");
		}
		// exclude the fully heterozygote vectors if the filter is set to true
		if (USE_FILTER_HETEROZYGOUS_FILTERING && genotypePattern.equals("ab/ab;ab/ab")) {
			throw new InvalidVCFLineException("Invalid VCF file: fully heterozygous variant", VCFLine); 
		}
		// exclude the 3/4 heterozygote vectors if the filter is set to true
		if (USE_FILTER_3QUATER_HETEROZYGOUS_FILTERING && 
				((genotypePattern.equals("ab+aa;ab/ab") || genotypePattern.equals("aa+ab;ab/ab") || genotypePattern.equals("ab/ab;aa/ab")))) {
			throw new InvalidVCFLineException("Invalid VCF file: 3/4 heterozygous variant", VCFLine); 
		}
	}


	/**
	 * Creates an instance of Variant with arbitrary informations for anything but genotype and inheritance state
	 * @param currentGenotype a genotype represented by a byte
	 */
	public Variant(int currentGenotype) {
		this.chromosome = "chr1";
		this.position = 1;
		this.referenceAllele = "C";
		this.alternatievAllele = "T";
		this.fatherAlleles = byteToAlleleType(currentGenotype, QuartetMember.FATHER);
		this.motherAlleles = byteToAlleleType(currentGenotype, QuartetMember.MOTHER);
		this.kid1Alleles = byteToAlleleType(currentGenotype, QuartetMember.KID1);
		this.kid2Alleles = byteToAlleleType(currentGenotype, QuartetMember.KID2);
		this.genotypePattern = computeGenotypePattern();
		this.quartetInheritanceStates = PatternToInheritanceStates.getInheritanceStates(genotypePattern);
		this.phasingQualityIndex = 0;
		//System.out.println(Integer.toBinaryString(currentGenotype) + "\t-->\t" + genotypePattern);
	}


	/**
	 * @param currentGenotype genotype represented as a byte with one bit per allele as follow:
	 * (paternal allele1, paternal allele2, maternal allele1, maternal allele2, kid1 allele1, kid1 allele2, kid2 allele1, kid2 allele2)
	 * @param member a {@link QuartetMember}
	 * @return the 2 alleles of the specified family member
	 */
	private AlleleType[] byteToAlleleType(int currentGenotype, QuartetMember member) {
		byte offsetAllele1 = 0;
		switch (member) {
		case FATHER:
			offsetAllele1 = 6;
			break;
		case MOTHER:
			offsetAllele1 = 4;
			break;
		case KID1:
			offsetAllele1 = 2;
			break;
		case KID2:
			offsetAllele1 = 0;
			break;
		}
		int offsetAllele2 = offsetAllele1 + 1;
		AlleleType allele1, allele2;
		if ((currentGenotype  & (1 << offsetAllele1)) == 0) {
			allele1 = AlleleType.REFERENCE_ALLELE;
		} else { 
			allele1 = AlleleType.ALTERNATIVE_ALLELE;
		}
		if ((currentGenotype  & (1 << offsetAllele2)) == 0) {
			allele2 = AlleleType.REFERENCE_ALLELE;
		} else { 
			allele2 = AlleleType.ALTERNATIVE_ALLELE;
		}
		AlleleType[] alleles = {allele1, allele2};
		return alleles;
	}


	/**
	 * @param genotypeFieldDesc the gentype field description from the VCF line (eg: GT:AD:DP:GQ:PL:PQ)
	 * @return the index of the PQ subfield (fields are separeted by ":") or -1 if not found
	 */
	private int getPhasingQualityFieldIndex(String genotypeFieldDesc) {
		String[] splitGenotypeFieldDesc = genotypeFieldDesc.split(":");
		for (int i = 0; i < splitGenotypeFieldDesc.length; i++) {
			if (splitGenotypeFieldDesc[i].trim().equals("PQ")) {
				return i;
			}
		}
		return -1;
	}


	/**
	 * @param genotypeField a genotype field from a VCF file
	 * @return true if the specified genotype field is phased
	 */
	private boolean isGenotypePhased(String genotypeField) {
		if (phasingQualityIndex != -1) {
			String[] splitGenotypeField = genotypeField.split(":");
			if (splitGenotypeField.length > phasingQualityIndex)  {
				double phasingQuality = Double.parseDouble(splitGenotypeField[phasingQualityIndex].trim());
				if (phasingQuality < PHASING_QUALITY_FILTER_THRESHOLD) {
					return false;
				}
			}
		}
		return (genotypeField.charAt(1) == '|');
	}


	/**
	 * 
	 * @param genotypeField format field following the GT:AD:DP:GQ:PL format
	 * @return the minimum PL score (ie the max probability that the genotype is not the one returned by the genotyper)
	 * @throws InvalidVCFLineException if there is more than 1 alternative
	 */
	private int genotypeFieldToPL(String genotypeField) throws InvalidVCFFieldException {
		String[] splitFormatField = genotypeField.split(":");
		// the following happens when we have a ./. variant
		if (splitFormatField.length < 5) {
			throw new InvalidVCFFieldException("Invalid VCF field: the genotype field has less than 5 subfield.", "Genotype Field", genotypeField);
		}
		String genotype = splitFormatField[0].trim();
		String[] plScores = splitFormatField[4].trim().split(",");
		int refRefScore = Integer.parseInt(plScores[0].trim());
		int refAltScore = Integer.parseInt(plScores[1].trim());
		int altAltScore = Integer.parseInt(plScores[2].trim());
		if ((genotype.equals("0/0")) || (genotype.equals("0|0"))) {
			return Math.min(refAltScore, altAltScore);
		}
		if ((genotype.equals("0/1")) || (genotype.equals("0|1")) || (genotype.equals("1|0"))) {
			return Math.min(refRefScore, altAltScore);
		}
		if ((genotype.equals("1/1")) || (genotype.equals("1|1"))) {
			return Math.min(refRefScore, refAltScore);
		}
		throw new InvalidVCFFieldException("Invalid VCF field.", "Genotype Field", genotypeField);
	}


	/**
	 * Throws an exception if the specified string is different from "PASS"
	 * @param filterField
	 * @throws FilteredVCFLineException
	 */
	private void filterFieldFiltering(String filterField) throws FilteredVCFLineException {
		switch (FILTER_FIELD_FILTERING) {
		case NINETY_NINE_POINT_NINE:
			if (!filterField.equalsIgnoreCase("PASS") && 
					!filterField.equalsIgnoreCase("TruthSensitivityTranche99.00to99.90") &&
					!filterField.equalsIgnoreCase("VQSRTrancheSNP99.90to100.00") &&
					!filterField.equalsIgnoreCase("VQSRTrancheINDEL99.90to100.00")) {
				throw new FilteredVCFLineException("Filter Field", filterField);
			}
			break;
		case PASS:
			if (!filterField.equalsIgnoreCase("PASS")) {
				throw new FilteredVCFLineException("Filter Field", filterField);
			}
			break;
		case NONE:
			// do nothing
			break;
		}
	}


	/**
	 * Hard filtering as recommended in the GATK best practice V3
	 * @param infoField info field of the VCF line
	 * @throws InvalidVCFLineException
	 */
	private void GATKFilter(String infoField) throws InvalidVCFFieldException, FilteredVCFLineException {
		// filter on the the QD field
		int QDIndex = infoField.indexOf("QD=");
		if (QDIndex == -1) {
			throw new InvalidVCFFieldException("Invalid VCF field: QD subfield not found", "Info Field", infoField);
		} else {
			String QDStr = infoField.substring(QDIndex + 3);
			int indexSemicolon = QDStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFFieldException("Invalid VCF field: QD subfield not found", "Info Field", infoField);
			} else {
				QDStr = QDStr.substring(0, indexSemicolon);
				double QD = Double.parseDouble(QDStr);
				//System.out.println("QD=" + QD);
				if (QD < 8.0) {
					throw new FilteredVCFLineException("QD", Double.toString(QD));
				}
			}
		}
		// filter on the the HRun field
		int HRunIndex = infoField.indexOf("HRun=");
		if (HRunIndex == -1) {
			throw new InvalidVCFFieldException("Invalid VCF field: HRun subfield not found", "Info Field", infoField);
		} else {
			String HRunStr = infoField.substring(HRunIndex + 5);
			int indexSemicolon = HRunStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFFieldException("Invalid HRun field: QD subfield not found", "Info Field", infoField);
			} else {
				HRunStr = HRunStr.substring(0, indexSemicolon);
				double HRun = Double.parseDouble(HRunStr);
				//System.out.println("HRun=" + HRun);
				if (HRun > 5) {
					throw new FilteredVCFLineException("HRun", Double.toString(HRun));
				}
			}			
		}
		// filter on the the FS field
		int FSIndex = infoField.indexOf("FS=");
		if (FSIndex == -1) {
			throw new InvalidVCFFieldException("Invalid VCF field: FS subfield not found", "Info Field", infoField);
		} else {
			String FSStr = infoField.substring(FSIndex + 3);
			int indexSemicolon = FSStr.indexOf(";");
			if (indexSemicolon == -1) {
				throw new InvalidVCFFieldException("Invalid VCF field: FS subfield not found", "Info Field", infoField);
			} else {
				FSStr = FSStr.substring(0, indexSemicolon);
				double FS = Double.parseDouble(FSStr);
				//System.out.println("FS=" + FS);
				if (FS > 200) {
					throw new FilteredVCFLineException("FS", Double.toString(FS));
				}			
			}
		}
	}


	/**
	 * @param genotypeInfo genotype information field of a VCF line
	 * @return the allele types of a sample
	 * @throws InvalidVCFFieldException
	 * @throws PartiallyCalledVariantException 
	 */
	private AlleleType[] stringToAlleleTypes(String genotypeInfo) throws InvalidVCFFieldException, PartiallyCalledVariantException {
		AlleleType[] resultAlleles = new AlleleType[2];
		String[] splitGenotypeInfo = genotypeInfo.split(":"); // the genotype info field is colon-separated
		String genotype = splitGenotypeInfo[0];	// the genotype is in the first field
		// the genotype is coded like follow 0/0, 0/1 (where 0 is the reference allele and 1 is the alternative allele)
		if (genotype.charAt(0) == '0') {
			resultAlleles[0] = AlleleType.REFERENCE_ALLELE;
		} else if (genotype.charAt(0) == '1') {
			resultAlleles[0] = AlleleType.ALTERNATIVE_ALLELE;
		} else if (genotype.charAt(0) == '.') {
			throw new PartiallyCalledVariantException(genotypeInfo);
		} else {
			throw new InvalidVCFFieldException("Invalid VCF field: the first allele value must be 0 or 1", "Genotype Field", genotypeInfo);
		}
		if(genotype.charAt(2) == '0') {
			resultAlleles[1] = AlleleType.REFERENCE_ALLELE;
		} else if (genotype.charAt(2) == '1') {
			resultAlleles[1] = AlleleType.ALTERNATIVE_ALLELE;
		} else if (genotype.charAt(0) == '.') {
			throw new PartiallyCalledVariantException(genotypeInfo);
		} else {
			throw new InvalidVCFFieldException("Invalid VCF field: the second allele value must be 0 or 1", "Genotype Field", genotypeInfo);
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
	 * The most frequent allele in the parents is denoted by �a�; in case of equal frequency, �a�
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
	 * @param allele an {@link AlleleType}
	 * @return the number of allele having the specified {@link AlleleType} in the quartet
	 */
	public final int getAlleleCount(AlleleType allele) {
		int alleleCount = 0;
		alleleCount = getAlleles(QuartetMember.FATHER)[0] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.FATHER)[1] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.MOTHER)[0] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.MOTHER)[1] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.KID1)[0] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.KID1)[1] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.KID2)[0] == allele ? alleleCount + 1 : alleleCount;
		alleleCount = getAlleles(QuartetMember.KID2)[1] == allele ? alleleCount + 1 : alleleCount;
		return alleleCount;
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
	public final String getAlternativeAllele() {
		return alternatievAllele;
	}


	/**
	 * @param member a quartet member
	 * @return the alleles of the specified quartet member
	 */
	public final AlleleType[] getAlleles(QuartetMember member) {
		switch (member) {
		case FATHER:
			return fatherAlleles;
		case MOTHER:
			return motherAlleles;
		case KID1:
			return kid1Alleles;
		case KID2:
			return kid2Alleles;
		default:
			return null;
		}
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
	public final QuartetInheritanceState[] getInheritanceStates() {
		return quartetInheritanceStates;
	}


	/**
	 * @param member a quartet member
	 * @return true if the specified quartet member is phased, false otherwise
	 */
	public final Boolean isPhased(QuartetMember member) {
		switch (member) {
		case FATHER:
			return isFatherPhased;
		case MOTHER:
			return isMotherPhased;
		case KID1:
			return isKid1Phased;
		case KID2:
			return isKid2Phased;
		default:
			return null;
		}
	}
	
	
	/**
	 * Set the phase of a specified quartet member
	 * @param member a quartet member
	 * @param isPhased true if the member's genotype is phased
	 */
	public void setPhase(QuartetMember member, boolean isPhased) {
		switch (member) {
		case FATHER:
			isFatherPhased = isPhased;
		case MOTHER:
			isMotherPhased = isPhased;
		case KID1:
			isKid1Phased = isPhased;
		case KID2:
			isKid2Phased = isPhased;
		}
	}


	/**
	 * Set the genotype of the specified quartet member
	 * @param member a quarte member
	 * @param firstAllele first allele of the member (paternal allele in children)
	 * @param secondAllele second allele of the member (maternal allele in children)
	 * @param isPhased true if the genotype is phased
	 */
	public final void setGenotype(QuartetMember member, AlleleType firstAllele, AlleleType secondAllele, boolean isPhased) {
		switch (member) {
		case FATHER:
			fatherAlleles[0] = firstAllele;
			fatherAlleles[1] = secondAllele;
			isFatherPhased = isPhased;
			break;
		case MOTHER:
			motherAlleles[0] = firstAllele;
			motherAlleles[1] = secondAllele;
			isMotherPhased = isPhased;
			break;
		case KID1:
			kid1Alleles[0] = firstAllele;
			kid1Alleles[1] = secondAllele;
			isKid1Phased = isPhased;
			break;
		case KID2:
			kid2Alleles[0] = firstAllele;
			kid2Alleles[1] = secondAllele;
			isKid2Phased = isPhased;
			break;
		}		
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
		variantString += quartetInheritanceStates[0];
		variantString += "\t";
		if (quartetInheritanceStates.length > 1) {
			variantString += quartetInheritanceStates[1];
		} else {
			variantString += "-";
		}
		return variantString;
	}


	/**
	 * Prints the variant on the bedgraph format
	 */
	public void printVariantBgrFormat() {
		System.out.println(this.getChromosome() + "\t" + this.getPosition() + "\t" + (this.getPosition() + 1) + "\t1");		
	}


	/**
	 * This methods prints in bgr format in the standard output the coordinates of the specified variant
	 * if the variant states correspond to the specified states
	 * @param states {@link QuartetInheritanceState}
	 */
	public void printVariantBgrFormat(QuartetInheritanceState... states) {
		int statesFoundCount = 0;		
		for (QuartetInheritanceState currentVariantState: this.getInheritanceStates()) {
			for (QuartetInheritanceState currentInputState: states) {
				if (currentVariantState == currentInputState) {
					statesFoundCount++;
				}
			}
		}
		if (statesFoundCount == states.length) {
			System.out.println(this.getChromosome() + "\t" + this.getPosition() + "\t" + this.getPosition() + "\t1");
		}
	}


	/**
	 * @param member a member of the family
	 * @return true if the variant is homozygous for the specified family member
	 */
	public boolean isHomozygous(QuartetMember member) {
		switch (member) {
		case FATHER:
			return fatherAlleles[0] == fatherAlleles[1];
		case MOTHER:
			return motherAlleles[0] == motherAlleles[1];
		case KID1:
			return kid1Alleles[0] == kid1Alleles[1];
		case KID2:
			return kid2Alleles[0] == kid2Alleles[1];
		default:
			return false;
		}
	}


	/**
	 * @param member a member of the family
	 * @return true if the variant is heterozygous for the specified family member
	 */
	public boolean isHeterozygous(QuartetMember member) {
		switch (member) {
		case FATHER:
			return fatherAlleles[0] != fatherAlleles[1];
		case MOTHER:
			return motherAlleles[0] != motherAlleles[1];
		case KID1:
			return kid1Alleles[0] != kid1Alleles[1];
		case KID2:
			return kid2Alleles[0] != kid2Alleles[1];
		default:
			return false;
		}
	}


	/**
	 * @param member a member of the family
	 * @return true if variant is a SNP for the specified member 
	 */
	public boolean isHomozygousReference(QuartetMember member) {
		switch (member) {
		case FATHER:
			return (fatherAlleles[0] == AlleleType.REFERENCE_ALLELE) && (fatherAlleles[1] == AlleleType.REFERENCE_ALLELE);
		case MOTHER:
			return (motherAlleles[0] == AlleleType.REFERENCE_ALLELE) && (motherAlleles[1] == AlleleType.REFERENCE_ALLELE);
		case KID1:
			return (kid1Alleles[0] == AlleleType.REFERENCE_ALLELE) && (kid1Alleles[1] == AlleleType.REFERENCE_ALLELE);
		case KID2:
			return (kid2Alleles[0] == AlleleType.REFERENCE_ALLELE) && (kid2Alleles[1] == AlleleType.REFERENCE_ALLELE);
		default:
			return false;
		}
	}


	/**
	 * @return true if the children are identical, false otherwise
	 */
	public boolean areChildrenIdentical() {
		// case where boths kids are phased
		if (isKid1Phased && isKid2Phased) {
			if ((kid1Alleles[0] == kid2Alleles[0]) && (kid1Alleles[1] == kid2Alleles[1])) {
				return true;
			}
		} else {
			// case where at least one kid is not phased
			if (((kid1Alleles[0] == kid2Alleles[0]) && (kid1Alleles[1] == kid2Alleles[1])) 
					|| ((kid1Alleles[0]) == kid2Alleles[1]) && (kid1Alleles[1] == kid2Alleles[0])) {
				return true;
			}
		}
		return false;
	}


	/**
	 * @return true if the Variant is an MIE, false otherwise
	 */
	public boolean isMIE() {
		for (QuartetInheritanceState currentState: quartetInheritanceStates) {
			if (currentState.equals(QuartetInheritanceState.MIE)) {
				return true;
			}
		}
		return false;
	}


	/**
	 * @return true if the Variant is not informative, false otherwise
	 */
	public boolean isNotInformative() {
		for (QuartetInheritanceState currentState: quartetInheritanceStates) {
			if (currentState.equals(QuartetInheritanceState.NOT_INFORMATIVE)) {
				return true;
			}
		}
		return false;
	}


	/**
	 * @return true if the variant is an indel
	 */
	public boolean isIndel() {
		String[] splitAlternativeAllele = alternatievAllele.split(",");
		int referenceLength = referenceAllele.length();
		boolean isIndel = false;
		int i = 0;
		while (!isIndel && i < splitAlternativeAllele.length) {
			isIndel = splitAlternativeAllele[i].trim().length() != referenceLength;
			i++;
		}
		return isIndel;
	}


	/**
	 * @return true if the variant is an insertion
	 */
	public boolean isInsertion() {
		String[] splitAlternativeAllele = alternatievAllele.split(",");
		int referenceLength = referenceAllele.length();
		boolean isInsertion = false;
		int i = 0;
		while (!isInsertion && i < splitAlternativeAllele.length) {
			isInsertion = splitAlternativeAllele[i].trim().length() > referenceLength;
			i++;
		}
		return isInsertion;
	}


	/**
	 * @return true if the variant is a deletion
	 */
	public boolean isDeletion() {
		String[] splitAlternativeAllele = alternatievAllele.split(",");
		int referenceLength = referenceAllele.length();
		boolean isDeletion = false;
		int i = 0;
		while (!isDeletion && i < splitAlternativeAllele.length) {
			isDeletion = splitAlternativeAllele[i].trim().length() < referenceLength;
			i++;
		}
		return isDeletion;

	}


	/**
	 * @return true if the variant is a SNP
	 */
	public boolean isSNP() {
		String[] splitAlternativeAllele = alternatievAllele.split(",");
		int referenceLength = referenceAllele.length();
		boolean isSNP = false;
		int i = 0;
		while (!isSNP && i < splitAlternativeAllele.length) {
			isSNP = splitAlternativeAllele[i].trim().length() == referenceLength;
			i++;
		}
		return isSNP;
	}


	/**
	 * @param vcfLine the vcf line that needs to be corrected
	 * @return the corrected paternal genotype
	 */
	public String getContaminationCorrectedGenotype(String vcfLine) {
		String[] splitLine = vcfLine.split("\t");
		String paternalGenotypeInfoField[] = splitLine[9].trim().split(":");
		String[] plScores = paternalGenotypeInfoField[4].trim().split(",");
		int refRefScore = Integer.parseInt(plScores[0].trim());
		int altAltScore = Integer.parseInt(plScores[2].trim());
		// if no alternative genotypes has a score of 100 or less we can't correct the paternal genotype
		if (Math.min(refRefScore, altAltScore) <= 80) {
			// case where the paternal genotype is 0|0 - this means that the maternal genotype
			// is 1|1 otherwise it wouldn't create a contamination
			if ((refRefScore < altAltScore) 
					&& ((getAlleles(QuartetMember.MOTHER)[0] == AlleleType.ALTERNATIVE_ALLELE)
							|| (getAlleles(QuartetMember.MOTHER)[1] == AlleleType.ALTERNATIVE_ALLELE))) {
				return "0/0";
			}
			// case where the paternal is 1|1
			if ((refRefScore > altAltScore) 
					&& ((getAlleles(QuartetMember.MOTHER)[0] == AlleleType.REFERENCE_ALLELE)
							|| (getAlleles(QuartetMember.MOTHER)[1] == AlleleType.REFERENCE_ALLELE))) {
				return "1/1";
			}
		}
		return null;
	}


	/**
	 * @param inheritanceState a {@link CrossTriosInheritanceState}
	 * @return true if the variant is a SCE for the specified {@link CrossTriosInheritanceState}
	 */
	public boolean isSCE(CrossTriosInheritanceState inheritanceState) {
		boolean isSCE = true;
		for (QuartetInheritanceState currentInheritanceState: getInheritanceStates()) {
			if (currentInheritanceState == QuartetInheritanceState.MIE) {
				return false;
			}
			if (currentInheritanceState == QuartetInheritanceState.NOT_INFORMATIVE) {
				return false;
			}
			if (currentInheritanceState.isCompatibleWith(inheritanceState)) {
				isSCE = false;
			}
		}
		return isSCE;
	}
}
