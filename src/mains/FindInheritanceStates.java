package mains;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import dataStructures.Variant;
import dataStructures.VariantListAnalyzer;

import exceptions.InvalidVCFLineException;


/**
 * Finds the inheritance states of a family quartet using the inheritance states as described in
 * "Analysis of Genetic Inheritance in a Family Quartet by Whole-Genome Sequencing", Roach et Al
 * @author Julien Lajugie
 */
public class FindInheritanceStates {

	// C:\Documents and Settings\Administrator\My Documents\GenPlay Library\Ritu_VCF\Ritu-corrected-ALL-LIBRARIES-SNP-chr1.raw.vcf

	/**
	 * Usage: java FindInheritanceStates.java -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java FindInheritanceStates.java -f <path to the file>");
			System.exit(-1);
		} else {
			try {
				findInheritancePatterns(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Finds the inheritance states of the family quartet
	 * @param VCFFilePath VCF files with the variants of the family quartet
	 * @throws IOException if the VCF file is not valid
	 */
	private static void findInheritancePatterns(String VCFFilePath) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFilePath));
			String line = null;
			List<Variant> variantList = new ArrayList<Variant>();			
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process variants with more than one alternative allele
						/*if (currentVariant.getAlternatievAllele().length() != 1) {
							System.err.println("Variants with more than 1 alternative allele are not handled:");
							System.err.println(line);
						}*/
						// we don't process variants that are not informative or a mandelien inheritance state
						/*if ((currentVariant.getInheritanceStates()[0] != InheritanceState.NOT_INFORMATIVE) &&
								(currentVariant.getInheritanceStates()[0] != InheritanceState.MIE)) {
							variantList.add(currentVariant);
						}*/
						variantList.add(currentVariant);
						//System.out.println(line);
					} catch (InvalidVCFLineException e) {
						// do nothing
					}					
				}
			}
			// analyze the variant list
			VariantListAnalyzer analyzer = new VariantListAnalyzer(variantList);
			// write in the standard output the result. Can be redirected in a file
			// write the header
			System.out.println(Variant.variantHeader() + "\tI\tM\tF\tNot-ID\tNot-Inf\tNotMend\tI Avg\tM Avg\tF Avg\tNot-ID Avg\tNot-Inf Avg\tNotMend Avg\tDominant");
			// write the result data			
			for (int i = 0; i < variantList.size(); i++) {
				System.out.print(variantList.get(i) + "\t");
				System.out.println(analyzer.getIdenticals()[i] + "\t" +
						analyzer.getMaternalIdenticals()[i] + "\t" +
						analyzer.getPaternalIdenticals()[i] + "\t" +
						analyzer.getNonIdenticals()[i] + "\t" +
						analyzer.getNonInfomatives()[i] + "\t" +
						analyzer.getNonMendelians()[i] + "\t" +
						analyzer.getIdenticalAvgs()[i] + "\t" +
						analyzer.getMaternalIdenticalAvgs()[i] + "\t" +
						analyzer.getPaternalIdenticalAvgs()[i] + "\t" +
						analyzer.getNonIdenticalAvgs()[i] + "\t" +
						analyzer.getNonInfomativeAvgs()[i] + "\t" +
						analyzer.getNonMendelianAvgs()[i] + "\t" +
						analyzer.getDominantInheritanceStates()[i]);
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
