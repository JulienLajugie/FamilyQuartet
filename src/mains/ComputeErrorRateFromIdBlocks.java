package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.QuartetInheritanceState;
import dataStructures.SegmentalDuplicationList;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;
import exceptions.VCFException;



/**
 * Computes the error rate by computing how many times the kids are different in identical blocks 
 * @author Julien Lajugie
 */
public class ComputeErrorRateFromIdBlocks {


	/**
	 * Usage: java ComputeErrorRateFromIdBlocks.java -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 * @param args -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ComputeErrorRateFromIdBlocks.java -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>");
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
				computeErrorRateFromIdBlocks(VCFFile, blockFile, segDupFile);
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
	 * Computes the error rate by computing how many times the kids are different in identical blocks
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @param segDupFile bed or bgr file containing the segmental duplication. Variants in these regions will be excluded.  Can be null
	 * @throws IOException if the VCF file is not valid
	 */
	private static void computeErrorRateFromIdBlocks(File VCFFile, File blockFile, File segDupFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		SegmentalDuplicationList segDupList = null;
		if (segDupFile != null) {
			segDupList = new SegmentalDuplicationList();
			segDupList.loadBedOrBgr(segDupFile);
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			int variantCount = 0;
			int errorCount = 0;
			int MIEErrorCount = 0;
			int SCEErrorCount = 0;
			long identicalBlocksLength = countIdenticalBlockLength(blockList);
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						if (segDupList == null || !segDupList.isInSegmentalDuplication(currentVariant)) {
							if (!currentVariant.isIndel()) {
								if ((blockList.getBlock(currentVariant) != null) && (blockList.getBlock(currentVariant).getBlockState().getMaternalTrioState() == TrioInheritanceState.IDENTICAL) 
										&& (blockList.getBlock(currentVariant).getBlockState().getPaternalTrioState() == TrioInheritanceState.IDENTICAL)) {
									variantCount++;
									if (!currentVariant.areChildrenIdentical()) {
										if (currentVariant.isMIE()) {
										currentVariant.printVariantBgrFormat();
										System.out.println(currentVariant.getGenotypePattern());
										}
										SCEErrorCount = blockList.getBlock(currentVariant).isSCE(currentVariant) ? SCEErrorCount + 1 : SCEErrorCount;
										MIEErrorCount = (currentVariant.getInheritanceStates()[0] == QuartetInheritanceState.MIE) ? MIEErrorCount + 1 : MIEErrorCount;
										errorCount++;
									}
								}
							}
						}
					} catch (VCFException e) {
						// do nothing
					}
				}
			}
			System.out.println("Variant Count=" + variantCount + ", Error Count=" + errorCount + ", Error Rate=" + (errorCount / (double) variantCount * 100d) + '%');
			System.out.println("MIE Error Count=" + MIEErrorCount + ", MIE Error %=" + (MIEErrorCount / (double) errorCount * 100d) + '%');
			System.out.println("SCE Error Count=" + SCEErrorCount + ", SCE Error %=" + (SCEErrorCount / (double) errorCount * 100d) + '%');
			System.out.println("Identical blocks total length=" + identicalBlocksLength + ", Identical Blocks Error Rate=" + (errorCount / (double) identicalBlocksLength));
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
	private static final long countIdenticalBlockLength(InheritanceStateBlockList<CrossTriosInheritanceState> blockList) {
		Collection<List<InheritanceStateBlock<CrossTriosInheritanceState>>> blockListCollection = blockList.getBlocks().values();
		long identicalBlockLength = 0;
		for (List<InheritanceStateBlock<CrossTriosInheritanceState>> currentList: blockListCollection) {
			for (InheritanceStateBlock<CrossTriosInheritanceState> currentBlock: currentList) {
				if ((currentBlock.getBlockState().getMaternalTrioState() == TrioInheritanceState.IDENTICAL) && 
						(currentBlock.getBlockState().getPaternalTrioState() == TrioInheritanceState.IDENTICAL)) {
					identicalBlockLength += (currentBlock.getStopPosition() - currentBlock.getStartPosition());
				}
			}
		}
		return identicalBlockLength;
	}
}
