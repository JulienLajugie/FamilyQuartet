package exceptions;

/**
 * Exception thrown when a value doesn't pass a filter
 * @author Julien Lajugie
 */
public class FilteredVCFLineException extends VCFException {
	
	private static final long serialVersionUID = 2408515088708705603L;	// generated ID
	private final String filterName;	// name of the filter that caused the exception
	private final String filteredValue;	// value that didn't pass the filter

	
	/**
	 * Creates an instance of {@link FilteredVCFLineException}
	 * @param filterName name of the filter that caused the exception
	 * @param filteredValue value that didn't pass the filter
	 */
	public FilteredVCFLineException(String filterName, String filteredValue) {
		super("A VCF line has been rejected by a filter");
		this.filterName = filterName;
		this.filteredValue = filteredValue;
	}

	
	/**
	 * @return the name of the filter that caused the exception
	 */
	public String getFilterName() {
		return filterName;
	}

	
	/**
	 * @return the value that didn't pass the filter
	 */
	public String getFilteredValue() {
		return filteredValue;
	}
}
