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
import exceptions.VCFException;


/**
 * Adds a MIE, SCE and RDF flags to the info field of the VCF file when it applies  
 * @author Julien Lajugie
 */
public class MarkSceMieRdfVariants {

	private static final String MIE_INFO_HEADER = "##INFO=<ID=MIE,Number=0,Type=Flag,Description=\"Mandelian Inheritance Error\">"; // header for the MIE flag
	private static final String SCE_INFO_HEADER = "##INFO=<ID=SCE,Number=0,Type=Flag,Description=\"State Consistency Error " +
			"(marked using inheritance state blocks from ";					// header for the SCE flag
	private static final String RDF_INFO_HEADER = "##INFO=<ID=RDF,Number=0,Type=Flag,Description=\"Variant rejected by the " +
			"Read Depth Filter (marked using segmental duplications from ";	// header for the RDF flag
	private static final String INFO_HEADER_END = ")\">"; // end of the SCE and RDF header
	private static final int	INFO_FIELD_INDEX = 7; // index of the info field in the vcf file

	/**
	 * Usage: java MarkSceMieRdfVariants -v <path to the VCF file> -b <path to the block file>  -s <path to the segmental duplication file>
	 * @param args -v <path to the VCF file> -b <path to the block file> -s <path to the segmental duplication file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java MarkSceMieRdfVariants.java -v <path to the VCF file> -b <path to the block file> -s <path to the segmental duplication file>");
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
				markSCEMIESegDup(VCFFile, blockFile, segDupFile);
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
		if (args.length != 6) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b")) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v") && !args[4].equals("-v")) {
			return false;
		}
		// case with no -s parameter 
		if (!args[0].equals("-s") && !args[2].equals("-s") && !args[4].equals("-s")) {
			return false;
		}
		return true;
	}


	/**
	 * Adds a MIE or SCE flags to the info field of the VCF file for the SCE and MIEvariants
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files in a bgr format
	 * @throws IOException if the VCF file is not valid
	 */
	private static void markSCEMIESegDup(File VCFFile, File blockFile, File segDupFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		SegmentalDuplicationList segDupList = new SegmentalDuplicationList();
		segDupList.loadBedOrBgr(segDupFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			boolean infoStarted = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					if (line.startsWith("##INFO")) {
						infoStarted = true;
					} else if (infoStarted) {
						System.out.println(MIE_INFO_HEADER);
						System.out.println(SCE_INFO_HEADER + blockFile.getPath() + INFO_HEADER_END);
						System.out.println(RDF_INFO_HEADER + blockFile.getPath() + INFO_HEADER_END);
						infoStarted= false;
					}
					System.out.println(line);
				} else {
					try {
						Variant currentVariant = new Variant(line);
						InheritanceStateBlock<CrossTriosInheritanceState> isBlock = blockList.getBlock(currentVariant);
						// a variant cannot be both MIE and SCE at the same time
						if (currentVariant.isMIE()) {
							line = markLineAs(line, "MIE");
						} else if ((isBlock != null) && (currentVariant.isSCE(isBlock.getBlockState()))) {
							line = markLineAs(line, "SCE");
						}
						if ((segDupList != null) && (segDupList.isInSegmentalDuplication(currentVariant))) {
							line = markLineAs(line, "RDF");
						}
					} catch (VCFException e) {
						// do nothing
					} finally {
						System.out.println(line);
					}
				}
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * @param line a vcf line
	 * @param flag a flag to add to the info field
	 * @return a copy of the vcf line where the specified flag was added to the info field 
	 */
	private static String markLineAs(String line, String flag) {
		String[] splitLine = line.split("\t");
		splitLine[INFO_FIELD_INDEX] += ";" + flag; 
		String newLine = splitLine[0];
		for (int i = 1; i < splitLine.length; i++) {
			newLine += "\t" + splitLine[i];
		}
		return newLine;
	}
}
