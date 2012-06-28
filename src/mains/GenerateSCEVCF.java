package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.QuartetInheritanceState;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;
import exceptions.VCFException;

/**
 * Generates a VCF containing only the SCE variants   
 * @author Julien Lajugie
 */
public class GenerateSCEVCF {


	/**
	 * Usage: java GenerateSCEVCF -v <path to the VCF file> -b <path to the block file>
	 * @param args -v <path to the VCF file> -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateSCEVCF -v <path to the VCF file> -b <path to the block file>");
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
	 * Generates a VCF containing only the SCE variants
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @throws IOException if the VCF file is not valid
	 */
	private static void generateBlockStats(File VCFFile, File blockFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);	
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			int i = 0 ;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process indels
						if (!currentVariant.isIndel()) {
						/*if (blockList.getBlock(currentVariant).isSCE(currentVariant)) {
							System.out.println(line);
						}*/
							String[] splitLine = line.split("\t"); // VCF fields are tab delimited
							// filter using the filter field
							if (splitLine[6].trim().equals("PASS") && i < 900000) {
								if (!currentVariant.isMIE() && ((blockList.getBlock(currentVariant) == null) || !currentVariant.isSCE(blockList.getBlock(currentVariant).getBlockState()))) {
									System.out.println(line);
									i++;
									//currentVariant.printVariantBgrFormat();
								}
							}
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
