package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.QuartetMember;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Adds the phase set information to a vcf file
 * @author Julien Lajugie
 */
public class AddPhaseSet {

	private static final String PS_FORMAT_HEADER = "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">"; // header for phase set flag of the format field
	
	
	/**
	 * Usage: java AddPhaseSet -v <path to the VCF file> -b <path to the block file>
	 * @param args -v <path to the VCF file> -b <path to the block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java AddPhaseSet.java -v <path to the VCF file> -b <path to the block file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File blockFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					} else if (args[i].equals("-b")) {
						blockFile = new File(args[i + 1]);
					}			
				}
				addPhaseSet(VCFFile, blockFile);
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
	 * Adds the phase set information to a vcf file
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile block files in a bgr format
	 * @throws IOException if the VCF file is not valid
	 */
	private static void addPhaseSet(File VCFFile, File blockFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			boolean formatStarted = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					if (line.startsWith("##FORMAT")) {
						formatStarted = true;
					} else if (formatStarted) {
						System.out.println(PS_FORMAT_HEADER);
						formatStarted= false;
					}
					System.out.println(line);
				} else {
					try {
						Variant currentVariant = new Variant(line);
						InheritanceStateBlock<CrossTriosInheritanceState> isBlock = blockList.getBlock(currentVariant);
						if (isBlock != null) {
							boolean isAtLeastOneMemberPhased = false;
							for (QuartetMember member: QuartetMember.values()) {
								if (currentVariant.isPhased(member)) {
									isAtLeastOneMemberPhased = true;
								}
							}
							if (isAtLeastOneMemberPhased) {
								line = addPhaseSet(line, isBlock, currentVariant);
							}
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
	 * @param line vcf line 
	 * @param isBlock block of the current variant
	 * @param variant variant of the vcf line
	 * @return the vcf line with the phase set subfield of the genotype info field set
	 */
	private static String addPhaseSet(String line, InheritanceStateBlock<CrossTriosInheritanceState> isBlock, Variant variant) {
		String[] splitLine = line.split("\t");
		splitLine[8] += ":PS";
		for (QuartetMember member: QuartetMember.values()) {
			int currentMemberIndex = getMemberInfoFieldIndex(member);
			if (!variant.isPhased(member)) {
				splitLine[currentMemberIndex] += ":.";
			} else {
				String phaseSet = String.valueOf(isBlock.getStartPosition());
				splitLine[currentMemberIndex] += ":" + phaseSet;
			}
		}	
		String newLine = splitLine[0];
		for (int i = 1; i < splitLine.length; i++) {
			newLine += "\t" + splitLine[i];
		}
		return newLine;
	}


	/**
	 * @param member a {@link QuartetMember}
	 * @return the index of the info field of the specified member in the vcf file
	 */
	private static int getMemberInfoFieldIndex(QuartetMember member) {
		switch (member) {
		case FATHER:
			return 9;
		case MOTHER:
			return 10;
		case KID1:
			return 11;
		case KID2:
			return 12;
		default :
			throw new IllegalArgumentException("Not a valid quartet member: " + member);
		}
	}
}
