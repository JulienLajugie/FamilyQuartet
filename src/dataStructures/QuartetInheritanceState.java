package dataStructures;


/**
 * This enumeration represent the different inheritance states of a quartet with two parents (founders) and two children which are:
 * haploidentical paternal,
 * haploidentical maternal,
 * identical,
 * non-identical,
 * not informative or
 * Mendelian inheritance error
 * @author Julien Lajugie
 */
public enum QuartetInheritanceState implements InheritanceState {

	/**
	 * Mendelian Inheritance Error
	 */
	MIE ("Mendelian Inheritance Error", -1),

	/**
	 * Partial block
	 */
	PARTIAL ("partial", 0),

	/**
	 * Non-identical state
	 */
	NON_IDENTICAL ("Non-Identical", 1),

	/**
	 * Haploidentical paternal state
	 */
	PATERNAL ("Haploidentical Paternal", 2),

	/**
	 * Haploidentical maternal state
	 */
	MATERNAL ("Haploidentical Maternal", 3),

	/**
	 * Identical state
	 */
	IDENTICAL ("Identical", 4),

	/**
	 * Not informative state
	 */
	NOT_INFORMATIVE ("Not Informative", 5);



	private final String 	name;	// name of the state 
	private final int 		score;	// score of the state


	/**
	 * Creates an instance of {@link QuartetInheritanceState}
	 * @param name name of the state
	 */
	private QuartetInheritanceState(String name, int score) {
		this.name = name;
		this.score = score;
	}


	/**
	 * @return the name of the state
	 */
	public String getName() {
		return name;
	}


	@Override
	public String toString() {
		return name;
	}


	/**
	 * @param score score representing an inheritance state as follow:
	 * -1 -> MIE
	 * 0 -> partial
	 * 1 -> Non-Identical
	 * 2 -> Paternal
	 * 3 -> Maternal
	 * 4 -> Identical
	 * 5 -> Not Informative
	 * @return the {@link QuartetInheritanceState} associated to the specified score
	 */
	public static QuartetInheritanceState valueOf(int score) {
		switch (score) {
		case -1:
			return QuartetInheritanceState.MIE;
		case 1:
			return QuartetInheritanceState.NON_IDENTICAL;
		case 2: 
			return QuartetInheritanceState.PATERNAL;
		case 3:
			return QuartetInheritanceState.MATERNAL;
		case 4: 
			return QuartetInheritanceState.IDENTICAL;
		case 5:
			return QuartetInheritanceState.NOT_INFORMATIVE;
		default:
			return null;			
		}
	}


	/**
	 * @param binaryState a binary representation of an {@link QuartetInheritanceState}
	 * @return the {@link QuartetInheritanceState} of the binary state
	 */
	public static QuartetInheritanceState valueOfBinaryState(String binaryState) {
		if (binaryState.equals("0000")) {
			return QuartetInheritanceState.IDENTICAL;
		}
		if (binaryState.equals("0001")) {
			return QuartetInheritanceState.PATERNAL;
		}
		if (binaryState.equals("0010")) {
			return QuartetInheritanceState.MATERNAL;
		}
		if (binaryState.equals("0011")) {
			return QuartetInheritanceState.NON_IDENTICAL;
		}
		return null;
	}


	/**
	 * @return the score value associated to the inheritance state 
	 */
	public int getScore() {
		return score;
	}


	/**
	 * @param blockState a {@link TrioInheritanceState}
	 * @return true if the {@link QuartetInheritanceState} is compatible with the specified {@link TrioInheritanceState}
	 */
	public boolean isCompatibleWith(CrossTriosInheritanceState blockState) {
		if ((blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN) && (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN)) {
			return true;
		}
		switch(this) {
		case IDENTICAL:
			if (((blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) || (blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN)) &&
					((blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) || (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN))) {
				return true;
			} else {
				return false;
			}
		case NON_IDENTICAL:
			if (((blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) || (blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN)) &&
					((blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) || (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN))) {
				return true;
			} else {
				return false;
			}			
		case PATERNAL:
			if (((blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) || (blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN)) &&
					((blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) || (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN))) {
				return true;
			} else {
				return false;
			}
		case MATERNAL:	
			if (((blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) || (blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN)) &&
					((blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) || (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN))) {
				return true;
			} else {
				return false;
			}
		case NOT_INFORMATIVE:
			if ((blockState.getPaternalTrioState() == TrioInheritanceState.NOT_INFORMATIVE) && (blockState.getMaternalTrioState() == TrioInheritanceState.NOT_INFORMATIVE)) {
				return true;
			} else {
				return false;
			}
		case MIE:
			if ((blockState.getPaternalTrioState() == TrioInheritanceState.MIE) || (blockState.getMaternalTrioState() == TrioInheritanceState.MIE)) {
				return true;
			} else {
				return false;
			}
		case PARTIAL:
			if ((blockState.getPaternalTrioState() == TrioInheritanceState.UNKNOWN) || (blockState.getMaternalTrioState() == TrioInheritanceState.UNKNOWN)) {
				return true;
			} else {
				return false;
			}
		}
		return false;
	}
}
