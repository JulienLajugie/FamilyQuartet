package dataStructures;

/**
 * This class repressents the differents 
 * @author Julien Lajugie
 */
public class CrossTriosInheritanceState implements InheritanceState {

	private final TrioInheritanceState paternalTrioState;	// inheritance state of the trio made of the father and two children
	private final TrioInheritanceState maternalTrioState;	// inheritance state of the trio made of the mother and two children


	/**
	 * Creates an instance of {@link CrossTriosInheritanceState}
	 * @param paternalTrioState inheritance state of the trio made of the father and two children
	 * @param maternalTrioState inheritance state of the trio made of the mother and two children
	 */
	public CrossTriosInheritanceState(TrioInheritanceState paternalTrioState, TrioInheritanceState maternalTrioState) {
		this.paternalTrioState = paternalTrioState;
		this.maternalTrioState = maternalTrioState;
	}


	/**
	 * @param score score score representing an inheritance state as follow:
	 * Digital State	Paternal State	Maternal State
	 * 0	Unknown	Unknown
	 * 1	Non-Identical	Unknown
	 * 2	Identical	Unknown
	 * 3	Unknown	Non-Identical
	 * 4	Non-Identical	Non-Identical
	 * 5	Identical	Non-Identical
	 * 6	Unknown	Identical
	 * 7	Non-Identical	Identical
	 * 8	Identical	Identical
	 * @return the {@link CrossTriosInheritanceState} associated to the specified score
	 */
	public static CrossTriosInheritanceState valueOf(int score) {
		int paternalScore = score % 3;
		int maternalScore = score / 3;
		TrioInheritanceState paternalTrioState = TrioInheritanceState.valueOf(paternalScore);
		TrioInheritanceState maternalTrioState = TrioInheritanceState.valueOf(maternalScore);
		return new CrossTriosInheritanceState(paternalTrioState, maternalTrioState);
	}
	

	@Override
	public String getName() {
		return paternalTrioState.getName() + " Paternal, " + maternalTrioState.getName() + " Maternal";
	}
	
	
	@Override
	public String toString() {
		return this.getName();
	}


	@Override
	public int getScore() {
		int score = paternalTrioState.getScore() + (maternalTrioState.getScore() * 3);
		return score;
	}
	
	
	/**
	 * @return the inheritance state of the trio made of the father and two children
	 */
	public final TrioInheritanceState getPaternalTrioState() {
		return paternalTrioState;
	}


	/**
	 * @return the inheritance state of the trio made of the mother and two children
	 */
	public final TrioInheritanceState getMaternalTrioState() {
		return maternalTrioState;
	}
	
	
	/**
	 * @param blockState a {@link QuartetInheritanceState}
	 * @return true if compatible with specified {@link QuartetInheritanceState}
	 */
	public boolean isCompatibleWith(QuartetInheritanceState blockState) {
		return blockState.isCompatibleWith(this);		
	}
}
