package dataStructures;

/**
 * This class represent a phased vector of a family quartet
 * @author Julien Lajugie
 *
 */
public class PhasedVector implements Comparable<PhasedVector> {
	
	private int 	position;	// position of the vector
	private short 	vector;		// genotype vector 
	
	
	/**
	 * Creates an instance of {@link PhasedVector}
	 * @param position position of the vector on a chromosome
	 * @param unphasedQuartetVector unphased vector represented by a string (ie: "abaaabab")
	 * @param phasedQuartetVector phased vector represented by a string (ie: "..aaabab")
	 */
	public PhasedVector(int position, String unphasedQuartetVector, String phasedQuartetVector) {
		this.position = position;
		this.vector = stringVector2IntVector(unphasedQuartetVector, phasedQuartetVector);
	}

	
	/**
	 * 
	 * @param familyQuartetVector
	 * @return
	 */
	private short stringVector2IntVector(String unphasedQuartetVector, String phasedQuartetVector) {
		short vector = 0;
		// father
		if ((phasedQuartetVector.charAt(0) != '.') && (phasedQuartetVector.charAt(1) != '.')) {
			// case sample phased
			vector = setBit(vector, 0, phasedQuartetVector.charAt(0) == 'a');
			vector = setBit(vector, 1, phasedQuartetVector.charAt(1) == 'a');
			vector = setBit(vector, 2, true);			
		} else {
			// case sample unphased
			vector = setBit(vector, 0, unphasedQuartetVector.charAt(0) == 'a');
			vector = setBit(vector, 1, unphasedQuartetVector.charAt(1) == 'a');
			vector = setBit(vector, 2, false);			
		}
		// mother
		if ((phasedQuartetVector.charAt(2) != '.') && (phasedQuartetVector.charAt(3) != '.')) {
			// case sample phased
			vector = setBit(vector, 3, phasedQuartetVector.charAt(2) == 'a');
			vector = setBit(vector, 4, phasedQuartetVector.charAt(3) == 'a');
			vector = setBit(vector, 5, true);			
		} else {
			// case sample unphased
			vector = setBit(vector, 3, unphasedQuartetVector.charAt(2) == 'a');
			vector = setBit(vector, 4, unphasedQuartetVector.charAt(3) == 'a');
			vector = setBit(vector, 5, false);			
		}
		// father
		if ((phasedQuartetVector.charAt(4) != '.') && (phasedQuartetVector.charAt(5) != '.')) {
			// case sample phased
			vector = setBit(vector, 6, phasedQuartetVector.charAt(4) == 'a');
			vector = setBit(vector, 7, phasedQuartetVector.charAt(5) == 'a');
			vector = setBit(vector, 8, true);			
		} else {
			// case sample unphased
			vector = setBit(vector, 6, unphasedQuartetVector.charAt(4) == 'a');
			vector = setBit(vector, 7, unphasedQuartetVector.charAt(5) == 'a');
			vector = setBit(vector, 8, false);			
		}
		// father
		if ((phasedQuartetVector.charAt(6) != '.') && (phasedQuartetVector.charAt(7) != '.')) {
			// case sample phased
			vector = setBit(vector, 9, phasedQuartetVector.charAt(6) == 'a');
			vector = setBit(vector, 10, phasedQuartetVector.charAt(7) == 'a');
			vector = setBit(vector, 11, true);			
		} else {
			// case sample unphased
			vector = setBit(vector, 9, unphasedQuartetVector.charAt(6) == 'a');
			vector = setBit(vector, 10, unphasedQuartetVector.charAt(7) == 'a');
			vector = setBit(vector, 11, false);			
		}
		return vector;
	}
	
	
	/**
	 * Sets the specified bit of a short primitive to the specified value
	 * @param vector short primitive containing 16 bits
	 * @param bitPosition position of the bit to set
	 * @param value value to set
	 */
	private short setBit(short vector, int bitPosition, boolean value) {
	  if (value) {
		  vector |= (1 << bitPosition);
	  } else {
		  vector &= ~(1 << bitPosition);
	  }
	  return vector;
	}
	
	
	/**
	 * @param vector short primitive containing 16 bits
	 * @param bitPosition position of the bit
	 * @return true if the bit at the specified position is set, false otherwise
	 */
	private boolean isSet(short vector, int bitPosition) {
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
		if (isSet(vector, 0)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(vector, 2)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(vector, 1)) {
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
		if (isSet(vector, 3)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(vector, 5)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(vector, 4)) {
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
		if (isSet(vector, 6)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(vector, 8)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(vector, 7)) {
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
		if (isSet(vector, 9)) {
			result = "0";
		} else {
			result = "1";
		}
		if (isSet(vector, 11)) {
			result += "|";
		} else {
			result += "/";
		}
		if (isSet(vector, 10)) {
			result += "0";
		} else {
			result += "1";
		}
		return result;	
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
}
