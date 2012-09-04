package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Computes statistics about phasing:
 * ie: number of phased and unphased variants, heterozygous variants per samples and number of quadruple heterozygous
 * @author Julien Lajugie
 */
public class ComputePhasingStats {


	/**
	 * Usage: java ComputePhasingStats.java -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java ComputePhasingStats.java -f <path to the file>");
			System.exit(-1);
		} else {
			File vcfFile = new File(args[1]);
			try {
				computePhasingStats(vcfFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Computes statistics about phasing
	 * @param vcfFile input vcf file
	 */
	private static void computePhasingStats(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			long fatherUnphased = 0;
			long fatherPhased = 0;
			long fatherHeteroPhased = 0;
			long motherUnphased = 0;
			long motherPhased = 0;
			long motherHeteroPhased = 0;
			long kid1Unphased = 0;
			long kid1Phased = 0;
			long kid1HeteroPhased = 0;
			long kid2Unphased = 0;
			long kid2Phased = 0;
			long kid2HeteroPhased = 0;
			long quadrupleHeterozygous = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						if (currentVariant.isSNP()) {
							if (currentVariant.isPhased(QuartetMember.FATHER)) {
								fatherPhased++;
								if (currentVariant.isHeterozygous(QuartetMember.FATHER)) {
									fatherHeteroPhased++;
								}
							} else {
								fatherUnphased++;
							}
							if (currentVariant.isPhased(QuartetMember.MOTHER)) {
								motherPhased++;
								if (currentVariant.isHeterozygous(QuartetMember.MOTHER)) {
									motherHeteroPhased++;
								}
							} else {
								motherUnphased++;
							}
							if (currentVariant.isPhased(QuartetMember.KID1)) {
								kid1Phased++;
								if (currentVariant.isHeterozygous(QuartetMember.KID1)) {
									kid1HeteroPhased++;
								}
							} else {
								kid1Unphased++;
							}
							if (currentVariant.isPhased(QuartetMember.KID2)) {
								kid2Phased++;
								if (currentVariant.isHeterozygous(QuartetMember.KID2)) {
									kid2HeteroPhased++;
								}
							} else {
								kid2Unphased++;
							}
							if (currentVariant.getGenotypePattern().equals("ab/ab;ab/ab")) {
								quadrupleHeterozygous++;
							}
						}
					} catch (VCFException e) {
						// do nothing
					}
				}
			}
			System.out.println("Sample\tunphased\tphased\thetero_phased");
			System.out.println("2_2\t" + fatherUnphased + "\t\t" + fatherPhased + "\t" + fatherHeteroPhased);
			System.out.println("2_5\t" + motherUnphased + "\t\t" + motherPhased + "\t" + motherHeteroPhased);
			System.out.println("3_2\t" + kid1Unphased + "\t\t" + kid1Phased + "\t" + kid1HeteroPhased);
			System.out.println("3_3\t" + kid2Unphased + "\t\t" + kid2Phased + "\t" + kid2HeteroPhased);
			System.out.println("Full Heterozygous: " + quadrupleHeterozygous);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
