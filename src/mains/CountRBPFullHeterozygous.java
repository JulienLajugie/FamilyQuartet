package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Counts the number of full heterozygous variants phased using the RBP algorithm 
 * @author Julien Lajugie
 */
public class CountRBPFullHeterozygous {
		
	/**
	 * Usage: java CountRBPFullHeterozygous.java -v <path RBP vcf file>
	 * @param args -v <path RBP vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java CountRBPFullHeterozygous.java -v <path RBP vcf file>");
			System.exit(-1);
		} else {
			File RBPFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-v")) {
					RBPFile = new File(args[i + 1]);
				}
			}
			try {
				countRBPFullHeterozygous(RBPFile);
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
	 * Counts the number of full heterozygous variants phased using the RBP algorithm
	 * @param RBPFile read backed phased VCF file
	 * @throws IOException
	 */
	private static void countRBPFullHeterozygous(File RBPFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(RBPFile));
			String line = null;
			int paternalPhasedVariantCount = 0;
			int maternalPhasedVariantCount = 0;
			int kid1PhasedVariantCount = 0;
			int kid2PhasedVariantCount = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						if (currentVariant.getGenotypePattern().equals("ab/ab;ab/ab")) {
							paternalPhasedVariantCount = currentVariant.isPhased(QuartetMember.FATHER) ? paternalPhasedVariantCount + 1 : paternalPhasedVariantCount;
							maternalPhasedVariantCount = currentVariant.isPhased(QuartetMember.MOTHER) ? maternalPhasedVariantCount + 1 : maternalPhasedVariantCount;
							kid1PhasedVariantCount = currentVariant.isPhased(QuartetMember.KID1) ? kid1PhasedVariantCount + 1 : kid1PhasedVariantCount;
							kid2PhasedVariantCount = currentVariant.isPhased(QuartetMember.KID2) ? kid2PhasedVariantCount + 1 : kid2PhasedVariantCount;
						}
					} catch (VCFException e) {
						// do nothing
					}					
				}
			}
			System.out.println("Paternal full heterozygous variants phased: " + paternalPhasedVariantCount);
			System.out.println("Maternal full heterozygous variants phased: " + maternalPhasedVariantCount);
			System.out.println("Kid1 full heterozygous variants phased: " + kid1PhasedVariantCount);
			System.out.println("Kid2 full heterozygous variants phased: " + kid2PhasedVariantCount);
			int totalPhasedVariants = paternalPhasedVariantCount + maternalPhasedVariantCount + kid1PhasedVariantCount + kid2PhasedVariantCount;
			System.out.println("Total full heterozygous variants phased: " + totalPhasedVariants);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}		
	}
}
