package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Counts the number of SNP variants for a each family member
 * @author Julien Lajugie
 */
public class CountSNPVariants {

	// C:\Documents and Settings\Administrator\My Documents\GenPlay Library\Ritu_VCF\Ritu-corrected-ALL-LIBRARIES-SNP-chr1.raw.vcf
	private static final long GENOME_LENGTH = 2897310462l; // length of hg19 genome (without N's)
	/**
	 * Usage: java CountSNPVariants.java -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java CountSNPVariants.java -f <path to the file>");
			System.exit(-1);
		} else {
			File vcfFile = new File(args[1]);
			try {
				countSNPVariants(vcfFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Counts the number of SNP variants for a each family member
	 * @param vcfFile input vcf file
	 */
	private static void countSNPVariants(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			long fatherSnpCount = 0;
			long motherSnpCount = 0;
			long kid1SnpCount = 0;
			long kid2SnpCount = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process indels variants
						if ((currentVariant.getReferenceAllele().length() == 1) && (currentVariant.getAlternatievAllele().length() == 1)) {
							if (currentVariant.isSNP(QuartetMember.FATHER)) {
								fatherSnpCount++;
							}
							if (currentVariant.isSNP(QuartetMember.MOTHER)) {
								motherSnpCount++;
							}
							if (currentVariant.isSNP(QuartetMember.KID1)) {
								kid1SnpCount++;
							}
							if (currentVariant.isSNP(QuartetMember.KID2)) {
								kid2SnpCount++;
							}
						}

					} catch (VCFException e) {
						// do nothing
					}	
				}
			}
			System.out.println("***FATHER*** Number of SNPs:" + fatherSnpCount + ", average distance between 2 SNPs:" + (GENOME_LENGTH / fatherSnpCount));
			System.out.println("***MOTHER*** Number of SNPs:" + motherSnpCount + ", average distance between 2 SNPs:" + (GENOME_LENGTH / motherSnpCount));
			System.out.println("***KID1*** Number of SNPs:" + kid1SnpCount + ", average distance between 2 SNPs:" + (GENOME_LENGTH / kid1SnpCount));
			System.out.println("***KID2*** Number of SNPs:" + kid2SnpCount + ", average distance between 2 SNPs:" + (GENOME_LENGTH / kid2SnpCount));
		} finally {
			if (reader != null) {
				reader.close();
			}
		}	
	}	
}
