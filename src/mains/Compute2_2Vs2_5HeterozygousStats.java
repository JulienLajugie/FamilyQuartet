package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Computes the number of positions where 2_2 is different from 2_5
 * @author Julien Lajugie
 */
public class Compute2_2Vs2_5HeterozygousStats {


	/**
	 * Usage: java Compute2_2Vs2_5HeterozygousStats.java -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java Compute2_2Vs2_5HeterozygousStats.java -f <path to the file>");
			System.exit(-1);
		} else {
			File vcfFile = new File(args[1]);
			try {
				compute2_2Vs2_5HeterozygousStats(vcfFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Computes the number of positions where 2_2 is different from 2_5
	 * @param vcfFile input vcf file
	 */
	private static void compute2_2Vs2_5HeterozygousStats(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			long SNPHomo = 0;
			long SNPHetero = 0;
			long InsHomo = 0;
			long InsHetero = 0;
			long DelHomo = 0;
			long DelHetero = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						if (currentVariant.isSNP()) {
							if ((currentVariant.getGenotypePattern().startsWith("aa+ab")) ||
									(currentVariant.getGenotypePattern().startsWith("ab+aa"))) {
								SNPHetero++;
							} else if (currentVariant.getGenotypePattern().startsWith("aa/bb")) {
								SNPHomo++;
							}
						}
						if (currentVariant.isInsertion()) {
							if ((currentVariant.getGenotypePattern().startsWith("aa+ab")) ||
									(currentVariant.getGenotypePattern().startsWith("ab+aa"))) {
								InsHetero++;
							} else if (currentVariant.getGenotypePattern().startsWith("aa/bb")) {
								InsHomo++;
							}
						}
						if (currentVariant.isDeletion()) {
							if ((currentVariant.getGenotypePattern().startsWith("aa+ab")) ||
									(currentVariant.getGenotypePattern().startsWith("ab+aa"))) {
								DelHetero++;
							} else if (currentVariant.getGenotypePattern().startsWith("aa/bb")) {
								DelHomo++;
							}
						}
					} catch (VCFException e) {
						// do nothing
					}
				}
			}
			System.out.println("***2_2 Vs 2_5***");
			System.out.println("***SNPs*** Number of heterozygous variants:" + SNPHetero);
			System.out.println("***SNPs*** Number of homozygous variants:" + SNPHomo);
			System.out.println("***Insertions*** Number of heterozygous variants:" + InsHetero);
			System.out.println("***Insertions*** Number of homozygous variants:" + InsHomo);
			System.out.println("***Deletions*** Number of heterozygous variants:" + DelHetero);
			System.out.println("***Deletions*** Number of homozygous variants:" + DelHomo);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
