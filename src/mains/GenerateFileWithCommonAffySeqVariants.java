package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.AffymetrixSNP;
import dataStructures.AffymetrixSNPList;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Creates a file containing the SNP that are common between the vcf from DNA-Seq and the SNP array from affymetrix.
 * The genotype in the result file is from the vcf DNA-Seq file
 * @author Julien Lajugie
 */
public class GenerateFileWithCommonAffySeqVariants {

	/**
	 * Usage: java GenerateFileWithCommonAffySeqVariants.java -v <path to the vcf file> -a <path to the affymetrix file>
	 * @param args -v <path to the vcf file> -a <path to the affymetrix file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateFileWithCommonAffySeqVariants.java -v <path to the vcf file> -a <path to the affymetrix file>");
			System.exit(-1);
		} else {
			File vcfFile = null;
			File affyFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-v")) {
					vcfFile = new File(args[i + 1]);
				}
				if (args[i].equals("-a")) {
					affyFile = new File(args[i + 1]);
				}
			}
			try {
				generateFileWithCommonAffySeqVariants(vcfFile, affyFile);
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
		// case with no -a parameter
		if (!args[0].equals("-a") && !args[2].equals("-a")) {
			return false;
		}
		return true;
	}


	/**
	 * Creates a file containing the SNP that are common between the vcf from DNA-Seq and the SNP array from affymetrix.
	 * The genotype in the result file is from the vcf DNA-Seq file
	 * @param vcfFile vcf file
	 * @param affyFile affymetrix SNP array text file
	 * @throws IOException
	 */
	private static void generateFileWithCommonAffySeqVariants(File vcfFile, File affyFile) throws IOException {
		// load affy file		
		AffymetrixSNPList affymetrixSNPList = new AffymetrixSNPList();
		affymetrixSNPList.loadAffymetrixFile(affyFile);

		//int VCFSNPFoundCount = 0;
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			System.out.println("#chromosome\tposition\tdbSNP\tVCF_genotype\tAffy_genotype");
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant variant = new Variant(line);
						// we don't want to work with indels
						if (!variant.isIndel()) {
							String chromosome = variant.getChromosome();
							int position = variant.getPosition();
							AffymetrixSNP affySNP = affymetrixSNPList.get(chromosome, position); 
							if (affySNP != null) {
								String[] splitLine = line.split("\t"); 
								String dbSNPRef = splitLine[2].trim();
								String vcfGenotype = vcfGenotypeToAffyGenotype(splitLine[9].trim());
								String affyGenotype = affySNP.getCallCode(); 
								System.out.println(chromosome + "\t" + position + "\t" + dbSNPRef + "\t" + vcfGenotype + "\t" + affyGenotype);
								//VCFSNPFoundCount++;
							}
						}
					} catch (VCFException e) {}					
				}
			}
			//System.out.println("Affy SNP Count: " + affymetrixSNPList.SNPCount());
			//System.out.println("VCF SNP Count: " + VCFSNPFoundCount);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * @param genotypeField genotype info field from a vcf file
	 * @return AA if the genotype is 0*0, AB if the genotype is 0*1 or 1*0, BB if the genotype is 1*1
	 */
	private static String vcfGenotypeToAffyGenotype(String genotypeField) {
		String genotype = genotypeField.split(":")[0].trim();
		if (genotype.charAt(0) == '0' && genotype.charAt(2) == '0') {
			return "AA";
		}
		if (genotype.charAt(0) == '1' && genotype.charAt(2) == '1') {
			return "BB";
		}
		return "AB";
	}
}
