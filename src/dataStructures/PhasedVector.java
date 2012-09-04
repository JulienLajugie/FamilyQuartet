package dataStructures;


/**
 * This class represent a phased vector of a family quartet
 * @author Julien Lajugie
 *
 */
public class PhasedVector implements Comparable<PhasedVector> {

	private int 	position;	// position of the vector
	private short 	vector;		// genotype vector as a short primitive containing 16 bits


	/**
	 * Creates an instance of {@link PhasedVector}
	 * @param position position of the vector on a chromosome
	 * @param variant variant to represent as a vector
	 */
	public PhasedVector(int position, Variant variant) {
		this.position = position;
		variant2Vector(variant);
	}
	
	
	/**
	 * Creates an instance of {@link PhasedVector}
	 * @param position position of the vector on a chromosome
	 * @param unphasedQuartetVector unphased vector represented by a string (eg: "abaaabab")
	 * @param phasedQuartetVector phased vector represented by a string (eg: "..aaabab")
	 */
	public PhasedVector(int position, String unphasedQuartetVector, String phasedQuartetVector) {
		this.position = position;
		stringVector2IntVector(unphasedQuartetVector, phasedQuartetVector);
	}


	/**
	 * Converts a variant to a to a vector
	 * @param variant variant to represent as a vector
	 */
	private void variant2Vector(Variant variant) {
		setBit(0, variant.getAlleles(QuartetMember.FATHER)[0] == AlleleType.REFERENCE_ALLELE);
		setBit(1, variant.getAlleles(QuartetMember.FATHER)[1] == AlleleType.REFERENCE_ALLELE);
		setBit(2, variant.isPhased(QuartetMember.FATHER));
		setBit(3, variant.getAlleles(QuartetMember.MOTHER)[0] == AlleleType.REFERENCE_ALLELE);
		setBit(4, variant.getAlleles(QuartetMember.MOTHER)[1] == AlleleType.REFERENCE_ALLELE);
		setBit(5, variant.isPhased(QuartetMember.MOTHER));
		setBit(6, variant.getAlleles(QuartetMember.KID1)[0] == AlleleType.REFERENCE_ALLELE);
		setBit(7, variant.getAlleles(QuartetMember.KID1)[1] == AlleleType.REFERENCE_ALLELE);
		setBit(8, variant.isPhased(QuartetMember.KID1));
		setBit(9, variant.getAlleles(QuartetMember.KID2)[0] == AlleleType.REFERENCE_ALLELE);
		setBit(10, variant.getAlleles(QuartetMember.KID2)[1] == AlleleType.REFERENCE_ALLELE);
		setBit(11, variant.isPhased(QuartetMember.KID2));
	}


	/**
	 * Converts a vector represented by a string (eg: "abaaabab") to a vector represented by a short
	 * @param unphasedQuartetVector unphased vector represented by a string (eg: "abaaabab")
	 * @param phasedQuartetVector phased vector represented by a string (eg: "..aaabab")
	 */
	private void stringVector2IntVector(String unphasedQuartetVector, String phasedQuartetVector) {
		// father
		if ((phasedQuartetVector.charAt(0) != '.') && (phasedQuartetVector.charAt(1) != '.')) {
			// case sample phased
			setBit(0, phasedQuartetVector.charAt(0) == 'a');
			setBit(1, phasedQuartetVector.charAt(1) == 'a');
			setBit(2, true);
		} else {
			// case sample unphased
			setBit(0, unphasedQuartetVector.charAt(0) == 'a');
			setBit(1, unphasedQuartetVector.charAt(1) == 'a');
			setBit(2, false);
		}
		// mother
		if ((phasedQuartetVector.charAt(2) != '.') && (phasedQuartetVector.charAt(3) != '.')) {
			// case sample phased
			setBit(3, phasedQuartetVector.charAt(2) == 'a');
			setBit(4, phasedQuartetVector.charAt(3) == 'a');
			setBit(5, true);
		} else {
			// case sample unphased
			setBit(3, unphasedQuartetVector.charAt(2) == 'a');
			setBit(4, unphasedQuartetVector.charAt(3) == 'a');
			setBit(5, false);
		}
		// father
		if ((phasedQuartetVector.charAt(4) != '.') && (phasedQuartetVector.charAt(5) != '.')) {
			// case sample phased
			setBit(6, phasedQuartetVector.charAt(4) == 'a');
			setBit(7, phasedQuartetVector.charAt(5) == 'a');
			setBit(8, true);
		} else {
			// case sample unphased
			setBit(6, unphasedQuartetVector.charAt(4) == 'a');
			setBit(7, unphasedQuartetVector.charAt(5) == 'a');
			setBit(8, false);
		}
		// father
		if ((phasedQuartetVector.charAt(6) != '.') && (phasedQuartetVector.charAt(7) != '.')) {
			// case sample phased
			setBit(9, phasedQuartetVector.charAt(6) == 'a');
			setBit(10, phasedQuartetVector.charAt(7) == 'a');
			setBit(11, true);
		} else {
			// case sample unphased
			setBit(9, unphasedQuartetVector.charAt(6) == 'a');
			setBit(10, unphasedQuartetVector.charAt(7) == 'a');
			setBit(11, false);
		}
	}


	/**
	 * Sets the specified bit of the vector
	 * @param bitPosition position of the bit to set
	 * @param value value to set
	 */
	private short setBit(int bitPosition, boolean value) {
		if (value) {
			vector |= (1 << bitPosition);
		} else {
			vector &= ~(1 << bitPosition);
		}
		return vector;
	}


	/**
	 * @param bitPosition position of the bit
	 * @return true if the bit at the specified position is set, false otherwise
	 */
	private boolean isSet(int bitPosition) {
		return (0x1 & (vector >> bitPosition)) == 1 ? true : false;
	}


	/**
	 * @return the position of the vector
	 */
	public int getPosition() {
		return position;
	}


	/**
	 * @return the father's phased haplotype
	 */
	public String getFatherGenotype() {
		String result = "";
		if (isSet(0)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(2)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(1)) {
			result += "0";
		} else {
			result += "1";
		}
		return result;
	}


	/**
	 * @return the mother's phased haplotype
	 */
	public String getMotherGenotype() {
		String result = "";
		if (isSet(3)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(5)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(4)) {
			result += "0";
		} else {
			result += "1";
		}
		return result;
	}


	/**
	 * @return the 1st kid's phased haplotype
	 */
	public String getKid1Genotype() {
		String result = "";
		if (isSet(6)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(8)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(7)) {
			result += "0";
		} else {
			result += "1";
		}
		return result;
	}


	/**
	 * @return the 2nd kid's phased haplotype
	 */
	public String getKid2Genotype() {
		String result = "";
		if (isSet(9)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(11)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(10)) {
			result += "0";
		} else {
			result += "1";
		}
		return result;	
	}


	/**
	 * @param quartetMember a {@link QuartetMember}
	 * @return the genotype of the specified family member
	 */
	public String getGenotype(QuartetMember quartetMember) {
		switch (quartetMember) {
		case FATHER:
			return getFatherGenotype();
		case MOTHER:
			return getMotherGenotype();
		case KID1:
			return getKid1Genotype();
		case KID2:
			return getKid2Genotype();
		default:
			return null;
		}
	}
	
	
	/**
	 * @param quartetMember  a {@link QuartetMember}
	 * @return true if the genotype of the specified member is phased. False otherwise
	 */
	public Boolean isPhased(QuartetMember quartetMember) {
		switch (quartetMember) {
		case FATHER:
			return isSet(2);
		case MOTHER:
			return isSet(5);
		case KID1:
			return isSet(8);
		case KID2:
			return isSet(11);
		default:
			return null;
		}
	}


	/**
	 * @param member a {@link QuartetMember}
	 * @return true if the specified member is heterozygous
	 */
	public Boolean isHeterozygous(QuartetMember member) {
		String genotype = getGenotype(member);
		if (genotype == null) {
			return null;
		} else {
			if (genotype.charAt(0) == genotype.charAt(2)) {
				return false;
			} else {
				return true;
			}
		}
	}
	
	
	/**
	 * Set the phasing of the specified family member to the specified value 
	 * @param quartetMember a quartet member
	 * @param isPhased true if the member is phased, false otherwise
	 */
	public void setPhasing(QuartetMember quartetMember, boolean isPhased) {
		int bitPosition = -1;
		switch (quartetMember) {
		case FATHER:
			bitPosition = 2;
			break;
		case MOTHER:
			bitPosition = 5;
			break;
		case KID1:
			bitPosition = 8;
			break;
		case KID2:
			bitPosition = 11;
			break;
		}
		if (bitPosition != -1) {
			setBit(bitPosition, isPhased);
		}
	}


	@Override
	public int compareTo(PhasedVector o) {
		if (this.getPosition() > o.getPosition()) {
			return 1;
		} else if (this.getPosition() < o.getPosition()) {
			return -1;
		} else {
			return 0;
		}
	}


	@Override
	public String toString() {
		return getFatherGenotype() + '+' + getMotherGenotype() + ';' + getKid1Genotype() + '+' + getKid2Genotype(); 
	}


	/**
	 * invert the phasing of the genotype of the specified member
	 * @param quartetMember
	 */
	public void invert(QuartetMember quartetMember) {
		int bitAllele1 = 0;
		int bitAllele2 = 0; 
		switch (quartetMember) {
		case FATHER:
			bitAllele1 = 0;
			bitAllele2 = 1;
			break;
		case MOTHER:
			bitAllele1 = 3;
			bitAllele2 = 4;			
			break;
		case KID1:
			bitAllele1 = 6;
			bitAllele2 = 7;
			break;
		case KID2:
			bitAllele1 = 9;
			bitAllele2 = 10;
			break;		
		}
		boolean bitTmp = isSet(bitAllele1);
		setBit(bitAllele1, isSet(bitAllele2));
		setBit(bitAllele2, bitTmp);		
	}
}
