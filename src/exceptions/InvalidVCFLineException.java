package exceptions;

/**
 * This exception is thrown when a line from a VCF file is not valid 
 * @author Julien Lajugie
 */
public class InvalidVCFLineException extends VCFException {

	private static final long serialVersionUID = 8080192983376331528L;
	private final String vcfLine;	// VCF line that caused the exception
	
	
	/**
	 * Creates an instance of an {@link InvalidVCFLineException}
	 * @param message detail message of the exception
	 * @param vcfLine VCF line that caused the exception
	 */
	public InvalidVCFLineException(String message, String vcfLine) {
		super(message);
		this.vcfLine = vcfLine;
	}


	/**
	 * @return the vcfLine that caused the exception
	 */
	public String getVCFLine() {
		return vcfLine;
	}
}
