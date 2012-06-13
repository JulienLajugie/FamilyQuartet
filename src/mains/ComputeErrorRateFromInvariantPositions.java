package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.AlleleType;
import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.SegmentalDuplicationList;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Computes the error rate by computing the errors positions that should be invariant aa/aa;aa/aa
 * @author Julien Lajugie
 */
public class ComputeErrorRateFromInvariantPositions {

	private static final long GENOME_LENGTH = 2897310462l; // length of hg19 genome (without N's)

	/**
	 * Usage: java ComputeErrorRateFromInvariantPositions.java -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 * @param args -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ComputeErrorRateFromInvariantPositions.java -v <path to the VCF file> -b <path to the block file> -s <segmental duplication file (optional)>");
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
				computeErrorRateFromInvariantPositions(VCFFile, blockFile, segDupFile);
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
	 * Computes the error rate by computing the errors positions that should be invariant aa/aa;aa/aa
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files from the ISCA software (Roach et Al)
	 * @param segDupFile bed or bgr file containing the segmental duplication. Variants in these regions will be excluded.  Can be null
	 * @throws IOException if the VCF file is not valid
	 */
	private static void computeErrorRateFromInvariantPositions(File VCFFile, File blockFile, File segDupFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		SegmentalDuplicationList segDupList = null;
		if (segDupFile != null) {
			segDupList = new SegmentalDuplicationList(segDupFile);
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			int variantCount = 0;
			int itemCount = 0;
			int errorCount = 0;
			int MIEErrorCount = 0;
			int SCEErrorCount = 0;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						if (segDupList == null || !segDupList.isInSegmentalDuplication(currentVariant)) {
							if (!currentVariant.isIndel()) {
								variantCount++;
								if ((currentVariant.getAlleleCount(AlleleType.REFERENCE_ALLELE) == 7) || (currentVariant.getAlleleCount(AlleleType.ALTERNATIVE_ALLELE) == 7)) {
									itemCount++;
									if (currentVariant.isMIE()) {
										MIEErrorCount++;
									} else if ((blockList.getBlock(currentVariant) != null) && (blockList.getBlock(currentVariant).isSCE(currentVariant))) {
										SCEErrorCount++;
									}
								}
							}
						}
					} catch (VCFException e) {
						// do nothing
					}				
				}
			}
			errorCount = MIEErrorCount + SCEErrorCount;
			long invariantPositionCount = GENOME_LENGTH - variantCount + errorCount;
			System.out.println("Variant Count=" + itemCount + ", Error Count=" + errorCount + ", Error Rate=" + (errorCount / (double) itemCount * 100d) + '%');
			System.out.println("MIE Error Count=" + MIEErrorCount + ", MIE Error %=" + (MIEErrorCount / (double) errorCount * 100d) + '%');
			System.out.println("SCE Error Count=" + SCEErrorCount + ", SCE Error %=" + (SCEErrorCount / (double) errorCount * 100d) + '%');
			System.out.println("Invariant Position Count=" + invariantPositionCount + ", Identical Blocks Error Rate per sample=" + (errorCount / (double) invariantPositionCount));
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
