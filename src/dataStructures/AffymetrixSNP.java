package dataStructures;

/**
 * A SNP from affymetrix SNP array
 * @author Julien Lajugie
 */
public class AffymetrixSNP implements Comparable<AffymetrixSNP> {

	private final String 	chromosome;		// chromosome of the SNP
	private final int 		position;		// position of the SNP
	private final String	callCode;		// affy genotype of the SNP
	
	/**
	 * creates an instance of {@link AffymetrixSNP}
	 * @param chromosome chromosome of the SNP
	 * @param position position of the SNP
	 * @param callCode affy genotype of the SNP
	 */
	public AffymetrixSNP(String chromosome, int position, String callCode) {
		this.chromosome = chromosome;
		this.position = position;
		this.callCode = callCode;		
	}
	
	
	/**
	 * @param affymetrixLine line from a affymetrix SNP array file
	 * @return an {@link AffymetrixSNP} if the it can be created from the line.  Null otherwise
	 */
	public static AffymetrixSNP createFromAffymetrixLine(String affymetrixLine) {
		String[] splitLine = affymetrixLine.split("\t");
		String callCode = splitLine[1].trim(); 
		if (splitLine.length < 9) {
			return null;
		}
		if (callCode.equalsIgnoreCase("nocall")) {
			return null;
		}
		String chromosome = "chr" + splitLine[7]; 
		int position = Integer.parseInt(splitLine[8]);
		return new AffymetrixSNP(chromosome, position, callCode);
	}

	
	/**
	 * @return the chromosome of the SNP
	 */
	public String getChromosome() {
		return chromosome;
	}
	

	/**
	 * @return the position of the SNP
	 */
	public int getPosition() {
		return position;
	}

	
	/**
	 * @return the affy genotype of the SNP
	 */
	public String getCallCode() {
		return callCode;
	}
	
	
	@Override
	public int compareTo(AffymetrixSNP o) {
		return Integer.valueOf(this.position).compareTo(Integer.valueOf(o.getPosition()));
	}
}
