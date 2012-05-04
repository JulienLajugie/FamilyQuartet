package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.QuartetInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;
import exceptions.VCFException;



/**
 * Computes the error rate by computing how many times the kids are different in identical blocks 
 * @author Julien Lajugie
 */
public class ComputeErrorRateFromIdBlocks {


	/**
	 * Usage: java ComputeErrorRateFromIdBlocks.java -v <path to the VCF file> -b <path to the block file>
	 * @param args -v <path to the VCF file> -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ComputeErrorRateFromIdBlocks -v <path to the VCF file> -b <path to the block file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = new File(args[0].equals("-v") ? args[1] : args[3]); 
				File blockFile = new File(args[0].equals("-b") ? args[1] : args[3]);
				computeErrorRateFromIdBlocks(VCFFile, blockFile);
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
	 * Computes the error rate by computing how many times the kids are different in identical blocks
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @throws IOException if the VCF file is not valid
	 */
	private static void computeErrorRateFromIdBlocks(File VCFFile, File blockFile) throws IOException {
		InheritanceStateBlockList<QuartetInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromQuartetBgrFile(blockFile);	
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			int variantCount = 0;
			int errorCount = 0;
			long identicalBlocksLength = countIdenticalBlockLength(blockList);
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
						if ((blockList.getBlock(currentVariant) != null) && (blockList.getBlock(currentVariant).getBlockState() == QuartetInheritanceState.IDENTICAL)) {
							variantCount++;
							if (currentVariant.getInheritanceStates()[0] == QuartetInheritanceState.MIE) {
								if (!currentVariant.areChildrenIdentical()) {
									System.out.println("MIE");
								}
							}
							if (blockList.getBlock(currentVariant).isSCE(currentVariant)) {
								if (!currentVariant.areChildrenIdentical()) {
									System.out.println("SCE");
								}
							}
							if (!currentVariant.areChildrenIdentical()) {
								errorCount++;
							}
						}
					} catch (VCFException e) {
						// do nothing
					}				
				}
			}
			System.out.println("Variant Count=" + variantCount + ", Error Count=" + errorCount + ", Error Rate=" + (errorCount / (double) variantCount * 100d) + '%');
			System.out.println("Identical blocks total length=" + identicalBlocksLength + ", Identical Blocks Error Rate=" + (errorCount / (double) identicalBlocksLength * 100d) + '%');
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
	
	
	/**
	 * @param blockList
	 * @return the sum of the length of the identical blocks
	 */
	private static final long countIdenticalBlockLength(InheritanceStateBlockList<QuartetInheritanceState> blockList) {
		Collection<List<InheritanceStateBlock<QuartetInheritanceState>>> blockListCollection = blockList.getBlocks().values();
		long identicalBlockLength = 0;
		for (List<InheritanceStateBlock<QuartetInheritanceState>> currentList: blockListCollection) {
			for (InheritanceStateBlock<QuartetInheritanceState> currentBlock: currentList) {
				if (currentBlock.getBlockState() == QuartetInheritanceState.IDENTICAL) {
					identicalBlockLength += (currentBlock.getStopPosition() - currentBlock.getStartPosition());
				}
			}
		}
		return identicalBlockLength;
	}
}
