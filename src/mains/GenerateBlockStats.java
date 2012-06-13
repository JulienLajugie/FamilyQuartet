package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.SegmentalDuplicationList;
import dataStructures.Variant;
import exceptions.PartiallyCalledVariantException;
import exceptions.VCFException;



/**
 * Generates statistics about the blocks of inheritance states (eg: count of MIE / SCE / not informative states per block)
 * "Analysis of Genetic Inheritance in a Family Quartet by Whole-Genome Sequencing", Roach et Al
 * @author Julien Lajugie
 */
public class GenerateBlockStats {


	/**
	 * Usage: java GenerateBlockStats -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 * @param args -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateBlockStats -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File blockFile = null;
				File segDupFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					} else if (args[i].equals("-b")) {
						blockFile = new File(args[i + 1]);
					} else if (args[i].equals("-s")) {			
						segDupFile = new File(args[i + 1]);
					}					
				}
				generateBlockStats(VCFFile, blockFile, segDupFile);
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
		if ((args.length != 4) && (args.length != 6)) {
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
	 * @param segDupFile bed or bgr file containing the segmental duplication. Variants in these regions will be excluded.  Can be null 
	 * @throws IOException if the VCF file is not valid
	 */
	private static void generateBlockStats(File VCFFile, File blockFile, File segDupFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		SegmentalDuplicationList segDupList = null;
		if (segDupFile != null) {
			segDupList = new SegmentalDuplicationList(segDupFile);
		}
		int variantCount = 0;
		int indelCount = 0;
		int partiallyCalledVariantCount = 0;
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
						if (segDupList == null || !segDupList.isInSegmentalDuplication(currentVariant)) {
							//if (currentVariant.isIndel()) {
								variantCount++;
								if (currentVariant.isIndel()) {
									indelCount++;
								}
								InheritanceStateBlock<CrossTriosInheritanceState> currentVariantBlock = blockList.getBlock(currentVariant);
								if (currentVariantBlock != null) {
									currentVariantBlock.analyzeVariant(currentVariant);
								}
							//}
						}
						//System.out.println(line);
					} catch (PartiallyCalledVariantException e) {
						// we still count partially called variants 
						variantCount++;
						partiallyCalledVariantCount++;
					} catch (VCFException e) {
						// do nothing
					}					
				}
			}	
			System.out.println("Variant #: " + variantCount + ", Partially Called Variant #: " + partiallyCalledVariantCount + ", SNP #: " + (variantCount - indelCount - partiallyCalledVariantCount) + ", Indel #: " + indelCount);
			blockList.printGenomeWideStatistics();
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
