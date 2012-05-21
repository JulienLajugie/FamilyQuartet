package dataStructures;

/**
 * This enumeration represents the different types of filter field from a VCF file
 * @author Julien Lajugie
 */
public enum VCFFilterField {

	/**
	 * NO Filter set
	 */
	NONE ("None"),
	
	/**
	 * Filter out the 99.90 to 100 tranche
	 */
	NINETY_NINE_POINT_NINE ("99.90"),
	
	/**
	 * Filter out the 99 to 100 tranche
	 */
	PASS ("PASS");
	
	private final String name;
	
	/**
	 * Creates an instance of {@link VCFFilterField}
	 * @param name name of the filter
	 */
	private VCFFilterField(String name) {
		this.name = name;
	}
	
	
	@Override
	public String toString() {
		return name;
	}
	
	
	/**
	 * @return the name of the filter
	 */
	public String getName() {
		return name;
	}
}
