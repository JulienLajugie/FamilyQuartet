package exceptions;

/**
 * Exception thrown when a variant is partially called (ie ./* genotype)
 * @author Julien Lajugie
 */
public class PartiallyCalledVariantException extends VCFException {

	private static final long serialVersionUID = 526554721820965692L; // generated ID
	private final String genotypeField;	// VCF line that caused the exception
	
	
	/**
	 * Creates an instance of an {@link PartiallyCalledVariantException}
	 * @param genotypeField genotype field that caused the exception
	 */
	public PartiallyCalledVariantException(String genotypeField) {
		super("The genotype of the variant is partially called");
		this.genotypeField = genotypeField;
	}


	/**
	 * @return the genotype field that caused the exception
	 */
	public String getGenotypeField() {
		return genotypeField;
	}
}
