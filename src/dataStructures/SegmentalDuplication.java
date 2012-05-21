package dataStructures;

/**
 * This class represents a segmental duplication define by a start and a stop position
 * @author Julien Lajugie
 */
public class SegmentalDuplication implements Comparable<SegmentalDuplication> {

	private final int 	startPosition;		// start position
	private final int 	stopPosition;		// stop position	
	
	
	/**
	 * Creates an instance of {@link SegmentalDuplication}
	 * @param startPosition start position of the segmental duplication
	 * @param stopPosition stop position of the segmental duplication
	 */
	public SegmentalDuplication(int startPosition, int stopPosition) {
		this.startPosition = startPosition;
		this.stopPosition = stopPosition;
	}

	
	/**
	 * @return the start position of the segmental duplication
	 */
	public int getStartPosition() {
		return startPosition;
	}
	
	
	/**
	 * @return the stop position of the segmental duplication
	 */
	public int getStopPosition() {
		return stopPosition;
	}


	@Override
	public int compareTo(SegmentalDuplication o) {
		int result = Integer.compare(this.getStartPosition(), o.getStartPosition());
		if (result != 0) {
			return result;
		} else {
			return Integer.compare(this.getStopPosition(), o.getStopPosition());		
		}
	}
}
