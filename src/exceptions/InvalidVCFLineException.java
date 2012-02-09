package exceptions;

public class InvalidVCFLineException extends Exception {

	private static final long serialVersionUID = 8080192983376331528L;

	public InvalidVCFLineException() {
		super("The Allele is not valid");
	}
}
