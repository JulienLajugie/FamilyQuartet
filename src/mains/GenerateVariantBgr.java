package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.InheritanceState;
import dataStructures.InheritanceStateBlockList;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;
import exceptions.VCFException;

/**
 * Prints the variants in a bgr format.  Only print the variant have the same state as the 
 * STATE_OF_VARIANTS_TO_PRINT constant. If this variable is set to null it will print the SCE variants. 
 * @author Julien Lajugie
 */
public class GenerateVariantBgr {


	private static final InheritanceState[] STATE_OF_VARIANTS_TO_PRINT = {InheritanceState.NON_IDENTICAL, InheritanceState.IDENTICAL};


	/**
	 * Usage: java GenerateVariantBgr -v <path to the VCF file> -b <path to the block file>
	 * @param args -v <path to the VCF file> -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateVariantBgr -v <path to the VCF file> -b <path to the block file>");
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
	 * Prints the variants in a bgr format.  Only print the variant have the same state as the 
	 * STATE_OF_VARIANTS_TO_PRINT constant. If this variable is set to null it will print the SCE variants. 
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
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process variants with more than one alternative allele or indels
						if ((currentVariant.getAlternatievAllele().length() != 1) || (currentVariant.getReferenceAllele().length() != 1)) {
							throw new InvalidVCFLineException("Invalid VCF line: indel or variant with more than one alt allele.", line);
						}
						if (STATE_OF_VARIANTS_TO_PRINT == null) {
							blockList.printSCEVariantBrgFormat(currentVariant);	
						} else {
							currentVariant.printVariantBgrFormat(STATE_OF_VARIANTS_TO_PRINT);
						}						
					} catch (VCFException e) {
						// do nothing
					}					
				}
			}	
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
