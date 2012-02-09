package dataStructures;


/**
 * This enumeration represent the different inheritance states which are:
 * haploidentical paternal,
 * haploidentical maternal,
 * identical,
 * non-identical,
 * not informative or
 * Mendelian inheritance error
 * @author Julien Lajugie
 */
public enum InheritanceState {
	
	/**
	 * Haploidentical paternal state
	 */
	PATERNAL ("Haploidentical Paternal"),
	
	/**
	 * Haploidentical maternal state
	 */
	MATERNAL ("Haploidentical Maternal"),
	
	/**
	 * Identical state
	 */
	IDENTICAL ("Identical"),
	
	/**
	 * Non-identical state
	 */
	NON_IDENTICAL ("Non-Identical"),
	
	/**
	 * Not informative state
	 */
	NOT_INFORMATIVE ("Not Informative"),
	
	/**
	 * Mendelian Inheritance Error
	 */
	MIE ("Mendelian Inheritance Error");
	
	
	private final String name;	// name of the state 
	
	
	/**
	 * Creates an instance of {@link InheritanceState}
	 * @param name name of the state
	 */
	private InheritanceState(String name) {
		this.name = name;
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
}
