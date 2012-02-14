package exceptions;

/**
 * This exception is thrown when a field from a line of a VCF file is not valid 
 * @author Julien Lajugie
 */
public class InvalidVCFFieldException extends VCFException {

	private static final long serialVersionUID = 8080192983376331528L;	// generated ID
	private final String vcfFieldName;	// name of the field that caused the exception
	private final String vcfFieldValue;	// value of the field that caused the exception
	
	
	/**
	 * Creates an instance of an {@link InvalidVCFFieldException}
	 * @param message detail message of the exception
	 * @param vcfFieldName name of the field that caused the exception
	 * @param vcfFieldValue value of the field that caused the exception
	 */
	public InvalidVCFFieldException(String message, String vcfFieldName, String vcfFieldValue) {
		super(message);
		this.vcfFieldName = vcfFieldName;
		this.vcfFieldValue = vcfFieldValue;
	}


	/**
	 * @return the name of the field that caused the exception
	 */
	public String getVCFFieldName() {
		return vcfFieldName;
	}


	/**
	 * @return the value of the field that caused the exception
	 */
	public String getVCFFieldValue() {
		return vcfFieldValue;
	}
}
