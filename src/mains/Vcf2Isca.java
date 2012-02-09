package mains;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * This class converts a VCF file to a haploscript input file 
 * @author Julien Lajugie <julien.lajugie@einstein.yu.edu>
 */
public class Vcf2Isca {

	// hard-coded headers of a haploscript file
	private static final String TAG_CREATOR = "#creator Julien Lajugie";
	private static final String TAG_SW_VERSION = "#sw_version 0.3-14-g8e4b69a (Mar 07 2011 13:28)";
	private static final String TAG_TIMESTAMP = "#timestamp " + new SimpleDateFormat("EEE MMM d HH:mm:ss yyyy").format(Calendar.getInstance().getTime());
	private static final String TAG_CONTACT = "#contact julien.lajugie@einstein.yu.edu";
	private static final String TAG_PEDIGREE_NAME = "#pedigree_name Quartet";
	private static final String TAG_PEDIGREE = "#pedigree 2_2 2_3 3_3:2_2,2_3 3_5:2_2,2_3";
	private static final String TAG_SEX = "#sex 2_2:m 2_3:f 3_3:f 3_5:f";
	private static final String TAG_PEDIGREE_VERSION = "#pedigree_version 20110501";
	private static final String TAG_PEDIGREE_DESCRIPTION = "#pedigree_description Bouhassira Lab Family Quartet Multigenerational Study";
	private static final String TAG_REFERENCE_GENOME = "#reference_genome hg19";
	private static final String TAG_NUMERICAL_BASE = "#numerical_base zero";
	private static final String TAG_ALLELE_REPRESENTATION = "#allele_representation refseq_abstraction";
	private static final String TAG_GENOTYPES = "#genotypes unphased";
	private static final String TAG_BLANK = "#";
	private static final String TAG_HEADER = "#header chromosome,start_position,reference,genotype_pattern,ind1_allele1_base,ind1_allele1_score,ind1_allele2_base,ind1_allele2_score,ind2_allele1_base,ind2_allele1_score,ind2_allele2_base,ind2_allele2_score,ind3_allele1_base,ind3_allele1_score,ind3_allele2_base,ind3_allele2_score,ind4_allele1_base,ind4_allele1_score,ind4_allele2_base,ind4_allele2_score";


	/**
	 * Usage: java Vcf2Isca -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if ((args == null) || 
				(args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java Vcf2Isca -f <path to the input vcf file>");
			System.exit(-1);
		} else {
			try {
				convert(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Convert a VCF file to a haploscript input file
	 * @param VCFFilePath
	 * @throws IOException
	 */
	private static void convert(String VCFFilePath) throws IOException {
		BufferedReader reader = null;
		try {
			// write the headers
			System.out.println(TAG_CREATOR);
			System.out.println(TAG_SW_VERSION);
			System.out.println(TAG_TIMESTAMP);
			System.out.println(TAG_CONTACT);
			System.out.println(TAG_PEDIGREE_NAME);
			System.out.println(TAG_PEDIGREE);
			System.out.println(TAG_SEX);
			System.out.println(TAG_PEDIGREE_VERSION);
			System.out.println(TAG_PEDIGREE_DESCRIPTION);
			System.out.println(TAG_REFERENCE_GENOME);
			System.out.println(TAG_NUMERICAL_BASE);
			System.out.println(TAG_ALLELE_REPRESENTATION);
			System.out.println(TAG_GENOTYPES);
			System.out.println(TAG_BLANK);
			System.out.println(TAG_HEADER);			
			// open the input file
			reader = new BufferedReader(new FileReader(VCFFilePath));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't want to take the comment lines into account
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					String position = splitLine[1].trim();
					String refAllele = splitLine[3].trim();
					String altAllele = splitLine[4].trim();
					String genotypePatern = retrieveGenotypePattern(splitLine);				
					String ind1Genotype = retrieveIndividualGenotype(splitLine[9].trim(), refAllele, altAllele);
					String ind2Genotype = retrieveIndividualGenotype(splitLine[10].trim(), refAllele, altAllele);
					String ind3Genotype = retrieveIndividualGenotype(splitLine[11].trim(), refAllele, altAllele);
					String ind4Genotype = retrieveIndividualGenotype(splitLine[12].trim(), refAllele, altAllele);				
					String outputLine = chromosome + "," + 
					position + "," +
					refAllele + "," +
					genotypePatern + "," +
					ind1Genotype + "," +
					ind2Genotype + "," +
					ind3Genotype + "," +
					ind4Genotype + ",";
					System.out.println(outputLine);	
				}
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * @param splitLine an array with all the field of the vcf file.
	 * The 10th to 13th fields contain the genotypes
	 * @return the genotype vector (eg: "aaababab")
	 */
	private static String retrieveGenotypePattern(String[] splitLine) {
		String ind1GenotypePattern = retrieveIndividualGenotypePattern(splitLine[9]);
		String ind2GenotypePattern = retrieveIndividualGenotypePattern(splitLine[10]);
		String ind3GenotypePattern = retrieveIndividualGenotypePattern(splitLine[11]);
		String ind4GenotypePattern = retrieveIndividualGenotypePattern(splitLine[12]);
		return ind1GenotypePattern + ind2GenotypePattern + ind3GenotypePattern + ind4GenotypePattern;
	}


	/**
	 * @param genotypeInfo genotype field of a vcf file 
	 * @return a string with the genotype pattern in the ISCA format ("aa", "ab", "ba", "bb")
	 */
	private static String retrieveIndividualGenotypePattern(String genotypeInfo) {
		String genotypePattern = genotypeInfo.split(":")[0];
		genotypePattern = genotypePattern.trim();
		genotypePattern = genotypePattern.replace('0', 'a');
		genotypePattern = genotypePattern.replace('1', 'b');
		genotypePattern = String.valueOf(genotypePattern.charAt(0)) + String.valueOf(genotypePattern.charAt(2));
		return genotypePattern;
	}


	/**
	 * @param genotypeInfo genotype field of the VCF file
	 * @param refAllele reference allele of the current variant
	 * @param altAllele alternative allele of the current variant
	 * @return a string with the 2 alleles and the 2 scores of one sample (eg: "A,33,T,33,")
	 */
	private static String retrieveIndividualGenotype(String genotypeInfo, String refAllele, String altAllele) {
		String[] splitGenotypeInfo = genotypeInfo.split(":"); // the genotype info field is colon-separated
		String score;
		if (splitGenotypeInfo.length > 3) {
			try {
				score = new Integer((int) Double.parseDouble(splitGenotypeInfo[3].trim())).toString();
			} catch (NumberFormatException e) {
				// case when the score is '.'
				score = "";
			}
		} else {
			score = "";
		}
		String genotype = splitGenotypeInfo[0].trim();
		String allele1;
		String allele2;
		if (genotype.charAt(0) == '0') {
			allele1 = refAllele;
		} else if (genotype.charAt(0) == '1') {
			allele1 = altAllele;
		} else {
			allele1 = "N";
		}
		if (genotype.charAt(2) == '0') {
			allele2 = refAllele;
		} else if (genotype.charAt(2) == '1') {
			allele2 = altAllele;
		} else {
			allele2 = "N";
		}
		return allele1 + "," + score + "," + allele2 + "," + score; 
	}
}
