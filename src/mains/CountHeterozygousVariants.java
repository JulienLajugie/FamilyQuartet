package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Counts the number of heterozygous variants for a each family member
 * @author Julien Lajugie
 */
public class CountHeterozygousVariants {

	// C:\Documents and Settings\Administrator\My Documents\GenPlay Library\Ritu_VCF\Ritu-corrected-ALL-LIBRARIES-SNP-chr1.raw.vcf
	private static final long GENOME_LENGTH = 2897310462l; // length of hg19 genome (without N's)
	/**
	 * Usage: java CountHeterozygousVariants.java -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java CountHeterozygousVariants.java -f <path to the file>");
			System.exit(-1);
		} else {
			File vcfFile = new File(args[1]);
			try {
				countHeterozygousVariants(vcfFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Counts the number of heterozygous variants for a each family member
	 * @param vcfFile input vcf file
	 */
	private static void countHeterozygousVariants(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			long fatherHeterozygousCount = 0;
			long motherHeterozygousCount = 0;
			long kid1HeterozygousCount = 0;
			long kid2HeterozygousCount = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process indels variants
						//if ((currentVariant.getReferenceAllele().length() == 1) && (currentVariant.getAlternatievAllele().length() == 1)) {
						
						// we want indels only
						//if ((currentVariant.getReferenceAllele().length() != 1) || (currentVariant.getAlternatievAllele().length() != 1)) {
						
							if (currentVariant.isHeterozygous(QuartetMember.FATHER)) {
								fatherHeterozygousCount++;
							}
							if (currentVariant.isHeterozygous(QuartetMember.MOTHER)) {
								motherHeterozygousCount++;
							}
							if (currentVariant.isHeterozygous(QuartetMember.KID1)) {
								kid1HeterozygousCount++;
							}
							if (currentVariant.isHeterozygous(QuartetMember.KID2)) {
								kid2HeterozygousCount++;
							}
						//}

					} catch (VCFException e) {
						// do nothing
					}	
				}
			}
			System.out.println("***FATHER*** Number of heterozygous variants:" + fatherHeterozygousCount + ", average distance between 2 variants:" + (GENOME_LENGTH / fatherHeterozygousCount));
			System.out.println("***MOTHER*** Number of heterozygous variants:" + motherHeterozygousCount + ", average distance between 2 variants:" + (GENOME_LENGTH / motherHeterozygousCount));
			System.out.println("***KID1*** Number of heterozygous variants:" + kid1HeterozygousCount + ", average distance between 2 variants:" + (GENOME_LENGTH / kid1HeterozygousCount));
			System.out.println("***KID2*** Number of heterozygous variants:" + kid2HeterozygousCount + ", average distance between 2 variants:" + (GENOME_LENGTH / kid2HeterozygousCount));
		} finally {
			if (reader != null) {
				reader.close();
			}
		}	
	}	
}
