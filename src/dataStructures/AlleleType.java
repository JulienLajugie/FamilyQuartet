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
	REFERENCE_ALLELE ("Reference Allele", 0),
	/**
	 * Alternative allele
	 */
	ALTERNATIVE_ALLELE ("Aleternative Allele", 1);

	private final String name;		// name of the allele type
	private final int	 intValue;	// integer value of an allele type

	/**
	 * Creates an instance of {@link AlleleType}
	 * @param name name of the allele type
	 * @param intValue integer value of an allele type
	 */
	private AlleleType(String name, int intValue) {
		this.name = name;
		this.intValue = intValue;
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

	
	/**
	 * @return the oposite allele (alternative if reference, reference if alternative)
	 */
	public AlleleType getOpositeAllele() {
		switch (this) {
		case REFERENCE_ALLELE:
			return ALTERNATIVE_ALLELE;
		case ALTERNATIVE_ALLELE:
			return REFERENCE_ALLELE;
		default:
			return null;
		}
	}


	/**
	 * @return the intValue associated to the {@link AlleleType}
	 */
	public int getIntValue() {
		return intValue;
	}
}
