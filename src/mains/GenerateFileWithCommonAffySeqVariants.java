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
			System.out.println("#chromosome\tposition\tdbSNP\tVCF_genotype\tAffy_genotype\tVCF_ref\tVCF_alt\tFilter\tmin_PL\tMIE\tSCE\tRDF");
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
								String vcfRef = variant.getReferenceAllele();
								String vcfAlt = variant.getAlternativeAllele();
								String filterField = splitLine[6];
								int minPL = getMinPl(splitLine);
								String infoField = splitLine[7];
								boolean isMIE = infoField.contains("MIE");
								boolean isSCE = infoField.contains("SCE");
								boolean isRDF = infoField.contains("RDF");
								System.out.println(chromosome + "\t" + position + "\t" + dbSNPRef + "\t" + vcfGenotype + "\t" + affyGenotype + "\t" + vcfRef + "\t" + vcfAlt + "\t" + filterField +"\t" + minPL + "\t" + isMIE + "\t" + isSCE + "\t" + isRDF);
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
	
	
	/**
	 * @param splitLine a line from a VCF file split by chromosome
	 * @return the minimum PL score of the samples of the VCF file
	 */
	private static int getMinPl(String[] splitLine) {
		int plScore = Math.min(genotypeFieldToPL(splitLine[9].trim()), genotypeFieldToPL(splitLine[10].trim()));
		plScore = Math.min(plScore, genotypeFieldToPL(splitLine[11].trim()));
		plScore = Math.min(plScore, genotypeFieldToPL(splitLine[12].trim()));
		return plScore;
	}
	
	
	/**
	 * @param genotypeField the genotype field of a sample
	 * @return the PL score of extracted from the specified genotype field
	 */
	private static int genotypeFieldToPL(String genotypeField) {
		String[] splitFormatField = genotypeField.split(":");
		String genotype = splitFormatField[0].trim();
		String[] plScores = splitFormatField[4].trim().split(",");
		int refRefScore = Integer.parseInt(plScores[0].trim());
		int refAltScore = Integer.parseInt(plScores[1].trim());
		int altAltScore = Integer.parseInt(plScores[2].trim());
		if ((genotype.equals("0/0")) || (genotype.equals("0|0"))) {
			return Math.min(refAltScore, altAltScore);
		}
		if ((genotype.equals("0/1")) || (genotype.equals("0|1")) || (genotype.equals("1|0"))) {
			return Math.min(refRefScore, altAltScore);
		}
		if ((genotype.equals("1/1")) || (genotype.equals("1|1"))) {
			return Math.min(refRefScore, refAltScore);
		}
		else return -1;
	}
}
