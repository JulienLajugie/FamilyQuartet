package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


/**
 * Generates a bedgraph file with the blocks from the read backed phasing VCF
 * @author Julien Lajugie
 */
public class GenerateBgrWithReadBackedBlocks {
		
	private static final int GENERATE_FATHER_BLOCKS = 0;
	@SuppressWarnings("unused")
	private static final int GENERATE_MOTHER_BLOCKS = 1;
	@SuppressWarnings("unused")
	private static final int GENERATE_KID1_BLOCKS = 2;
	@SuppressWarnings("unused")
	private static final int GENERATE_KID2_BLOCKS = 3;
	private static final int SELECTED_SAMPLE = GENERATE_FATHER_BLOCKS;

	/**
	 * Usage: java GenerateBgrWithReadBackedBlocks -v <path to the vcf file>
	 * @param args -v <path to the vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateBgrWithReadBackedBlocks -v <path to the vcf file>");
			System.exit(-1);
		} else {
			File vcfFile = new File(args[1]);
			try {
				generateBgrWithReadBackedBlocks(vcfFile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * @param args parameters from the main function
	 * @return true if the parameters are valid
	 */
	private static boolean areParametersValid(String[] args) {
		if (args == null) {
			return false;
		}
		if (args.length != 2) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v")) {
			return false;
		}
		return true;
	}


	/**
	 * Generates a bedgraph file with the blocks from the read backed phasing VCF
	 * @param vcfFile
	 */
	private static void generateBgrWithReadBackedBlocks(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			// loop until eof
			int startPositionTmp = -1;
			int startPosition = 0;
			int stopPosition = -1;
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					String genotype = splitLine[SELECTED_SAMPLE + 9].trim().split(":")[0];
					if (genotype.charAt(1) == '|') {
						stopPosition = position;
						if (startPositionTmp != -1) {
							startPosition = startPositionTmp;
							startPositionTmp = -1;
						}
					} else if (genotype.charAt(1) == '/') {
						startPositionTmp = position;
						if (stopPosition != -1) {
							System.out.println(chromosome + '\t' + startPosition + '\t' + stopPosition +"\t1");
							stopPosition = -1;
						}
					}
				}
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
