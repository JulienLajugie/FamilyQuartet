package dataStructures;

/**
 * Interface that needs to be implemented by the different type of inheritance states
 * @author Julien Lajugie
 */
public interface InheritanceState {
	
	/**
	 * @return the name of the state
	 */
	public String getName();
	
	
	/**
	 * @return the score associated to the state
	 */
	public int getScore();
}
