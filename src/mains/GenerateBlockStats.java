package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.InheritanceStateBlockList;
import dataStructures.Variant;

import exceptions.InvalidVCFLineException;
import exceptions.VCFException;



/**
 * Generates statistics about the blocks of inheritance states (eg: count of MIE / SCE / not informative states per block)
 * "Analysis of Genetic Inheritance in a Family Quartet by Whole-Genome Sequencing", Roach et Al
 * @author Julien Lajugie
 */
public class GenerateBlockStats {


	/**
	 * Usage: java GenerateBlockStats -v <path to the VCF file> -b <path to the block file>
	 * @param args -v <path to the VCF file> -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateBlockStats -v <path to the VCF file> -b <path to the block file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = new File(args[0].equals("-v") ? args[1] : args[3]); 
				File blockFile = new File(args[0].equals("-b") ? args[1] : args[3]);
				generateBlockStats(VCFFile, blockFile);
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
		if ((args[0].equals("-v") && args[2].equals("-b")) 
				|| (args[0].equals("-b") && args[2].equals("-v"))) {
			return true;
		} else {
			return false;
		}
	}


	/**
	 * Generates statistics about the blocks of inheritance states (eg: count of MIE / SCE / not informative states per block)
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @throws IOException if the VCF file is not valid
	 */
	private static void generateBlockStats(File VCFFile, File blockFile) throws IOException {
		InheritanceStateBlockList blockList = new InheritanceStateBlockList();
		blockList.loadFromISCAFile(blockFile);		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			//List<Variant> variantList = new ArrayList<Variant>();			
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process variants with more than one alternative allele or indels
						if ((currentVariant.getAlternatievAllele().length() != 1) || (currentVariant.getReferenceAllele().length() != 1)) {
							throw new InvalidVCFLineException("Indel or Variant with more than one alt allele", line);
						}
						blockList.analyzeVariant(currentVariant);
						//System.out.println(line);
					} catch (VCFException e) {
						// do nothing
					}					
				}
			}	
			blockList.printStatistics();
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
