package dataStructures;

/**
 * Enum representing the different type of alleles.
 * An allele can be a reference allele or a alternative allele
 * @author Julien Lajugie
 */
public enum AlleleType {
	
	/**
	 * Reference allele
	 */
	REFERENCE_ALLELE ("Reference Allele"),
	/**
	 * Alternative allele
	 */
	ALTERNATIVE_ALLELE ("Aleternative Allele");
	
	private final String name;	// name of the allele type
	
	
	/**
	 * Creates an instance of {@link AlleleType}
	 * @param name name of the allele type
	 */
	private AlleleType(String name) {
		this.name = name;
	}


	/**
	 * @return the name of the allele type
	 */
	public String getName() {
		return name;
	}
	
	
	@Override
	public String toString() {
		return name;
	}
}
