package dataStructures;


/**
 * This enumeration represent the different inheritance states of a trio with one parent (founder) and two children which are:
 * unknown,
 * Mendelian inheritance error,
 * identical,
 * non-identical,
 * not informative
 * @author Julien Lajugie
 */
public enum TrioInheritanceState implements InheritanceState {

	/**
	 * Mendelian Inheritance Error
	 */
	MIE ("Mendelian Inheritance Error", -1),
	
	/**
	 * Partial block
	 */
	UNKNOWN ("partial", 0),
	
	/**
	 * Non-identical state
	 */
	NON_IDENTICAL ("Non-Identical", 1),
	
	/**
	 * Identical state
	 */
	IDENTICAL ("Identical", 2),
	
	/**
	 * Not informative state
	 */
	NOT_INFORMATIVE ("Not Informative", 10);

		
	private final String 	name;	// name of the state
	private final int		score;	// score of the state
	
	
	/**
	 * Creates an instance of {@link TrioInheritanceState}
	 * @param name name of the state
	 * @score 
	 */
	private TrioInheritanceState(String name, int score) {
		this.name = name;
		this.score = score;
	}


	/**
	 * @return the name of the state
	 */
	@Override	
	public String getName() {
		return name;
	}
	
	
	@Override
	public String toString() {
		return name;
	}
	
	
	/**
	 * @param score score representing an inheritance state as follow:
	 * -2 -> not informative
	 * -1 -> MIE
	 * 0 -> Unknown
	 * 1 -> Non-identical
	 * 2 -> Identical
	 * @return the {@link TrioInheritanceState} associated to the specified score
	 */
	public static TrioInheritanceState valueOf(int score) {
		switch (score) {
		case -2:
			return TrioInheritanceState.NOT_INFORMATIVE;
		case -1:
			return TrioInheritanceState.MIE;
		case 0:
			return TrioInheritanceState.UNKNOWN;
		case 1:
			return TrioInheritanceState.NON_IDENTICAL;
		case 2: 
			return TrioInheritanceState.IDENTICAL;
		default:
			return null;			
		}
	}


	@Override
	public int getScore() {
		return score;
	}
}
