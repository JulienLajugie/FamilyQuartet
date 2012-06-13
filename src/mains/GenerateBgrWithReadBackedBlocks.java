package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Generates a bedgraph file with the blocks from the read backed phasing VCF
 * @author Julien Lajugie
 */
public class GenerateBgrWithReadBackedBlocks {

	/**
	 * Usage: java GenerateBgrWithReadBackedBlocks.java -v <path to the vcf file> -m <quartet member (FATHER, MOTHER, KID1 or KID2)>
	 * @param args -v <path to the vcf file> -m <quartet member (FATHER, MOTHER, KID1 or KID2)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("java GenerateBgrWithReadBackedBlocks.java -v <path to the vcf file> -m <quartet member (FATHER, MOTHER, KID1 or KID2)>");
			System.exit(-1);
		} else {
			QuartetMember member = null;
			File vcfFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-m")) {
					member = QuartetMember.valueOf(args[i + 1]);
					if (member == null) {
						System.out.println("Usage: java GenerateBgrWithReadBackedBlocks.java -v <path to the vcf file> -m <quartet member (FATHER, MOTHER, KID1 or KID2)>");
						System.exit(-1);						
					}
				}
				if (args[i].equals("-v")) {
					vcfFile = new File(args[i + 1]);
				}
			}
			try {
				generateBgrWithReadBackedBlocks(vcfFile, member);
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
		if (args.length != 4) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v")) {
			return false;
		}
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m")) {
			return false;
		}
		return true;
	}


	/**
	 * Generates a bedgraph file with the blocks from the read backed phasing VCF
	 * @param vcfFile input phased VCF file
	 * @param member member of the quartet 
	 */
	private static void generateBgrWithReadBackedBlocks(File vcfFile, QuartetMember member) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			// loop until eof
			int startPosition = 0;
			int stopPosition = 0;
			String chromosome = null;
			boolean isBlockFirstVariant = true;			
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					Variant currentVariant;
					try {
						currentVariant = new Variant(line);						
						// case where the variant is not phased or is on a new chromosome (meaning that the previous block ended)
						if ((!currentVariant.isPhased(member)) || (!currentVariant.getChromosome().equals(chromosome))) {
							// case where the previous block is not empty
							if (!isBlockFirstVariant) {
								System.out.println(chromosome + '\t' + startPosition + '\t' + (stopPosition + 1) +"\t1");
								isBlockFirstVariant = true;
							}
							chromosome = currentVariant.getChromosome();
							startPosition = currentVariant.getPosition();
						}
						// case where the variant is phased with the previsous one
						if ((currentVariant.isPhased(member)) && (currentVariant.isHeterozygous(member))) {
							isBlockFirstVariant = false;
							stopPosition = currentVariant.getPosition();
						} 
					} catch (VCFException e) {
						// do nothing
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
