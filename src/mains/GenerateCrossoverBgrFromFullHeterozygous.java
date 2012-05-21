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
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;
import exceptions.FilteredVCFLineException;
import exceptions.InvalidVCFFieldException;
import exceptions.InvalidVCFLineException;
import exceptions.PartiallyCalledVariantException;



/**
 * Generates a bedgraph showing the identical and non identical variants for a founder of the family quartet
 * using the block information as well as the full heterozygous variants
 * @author Julien Lajugie
 */
public class GenerateCrossoverBgrFromFullHeterozygous {

	/**
	 * Usage: java GenerateCrossoverBgrFromFullHeterozygous.java -v <path to the merge vcf file> -b <block bgr file> -m <FATHER or MOTHER>
	 * @param args -v <path to the merge vcf file> -b <block bgr file> -m <FATHER or MOTHER>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateCrossoverBgrFromFullHeterozygous.java -v <path to the merge vcf file> -b <block bgr file> -m <FATHER or MOTHER>");
			System.exit(-1);
		} else {
			QuartetMember member = null;
			File vcfFile = null;
			File blockFile = null;
			for (int i = 0; i <= 4; i += 2) {
				if (args[i].equals("-m")) {
					member = QuartetMember.valueOf(args[i + 1]);
					if ((member == null) || ((member != QuartetMember.FATHER) && (member != QuartetMember.MOTHER))) {
						System.out.println("Usage: java GenerateCrossoverBgrFromFullHeterozygous.java -v <path to the merge vcf file> -b <block bgr file> -m <FATHER or MOTHER>");
						System.exit(-1);						
					}
				}
				if (args[i].equals("-v")) {
					vcfFile = new File(args[i + 1]);
				}
				if (args[i].equals("-b")) {
					blockFile = new File(args[i + 1]);
				}
			}
			try {
				generateCrossoverBgrFromFullHeterozygous(vcfFile, blockFile, member);
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
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m") && !args[4].equals("-m")) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v") && !args[4].equals("-v")) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b")) {
			return false;
		}		
		return true;
	}


	/**
	 * Generate a bedgraph showing the identical and non identical variants for a founder of the family quartet
	 * @param vcfFile input vcf file with the variants for the quartet
	 * @param blockFile bgr file with the inheritance state blocks of the family
	 * @param founderMember a founder member of the family quartet
	 * @throws IOException
	 */
	private static void generateCrossoverBgrFromFullHeterozygous(File vcfFile, File blockFile, QuartetMember founderMember) throws IOException {
		// load the block file
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		int founderAllele = (founderMember == QuartetMember.MOTHER ? 0 : 1); // the 0 is the index of the paternal allele, 1 maternal
		BufferedReader reader = null;		
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant variant = new Variant(line);
						if (variant.getGenotypePattern().equals("ab/ab;ab/ab")) {
							int score;
							if (founderMember == QuartetMember.FATHER) {
								InheritanceStateBlock<CrossTriosInheritanceState> block = blockList.getBlock(variant);
								if ((block != null) && (block.getBlockState().getMaternalTrioState() == TrioInheritanceState.IDENTICAL)) {
									score = 1;
								} else {
									score = -1;
								}
							} else {
								// case where we generate the bgr for the mother 
								InheritanceStateBlock<CrossTriosInheritanceState> block = blockList.getBlock(variant);
								if ((block != null) && (block.getBlockState().getPaternalTrioState() == TrioInheritanceState.IDENTICAL)) {
									score = 1;
								} else {
									score = -1;
								}								
							}					
							System.out.println(variant.getChromosome() + '\t' + variant.getPosition() + '\t' + (variant.getPosition() + 1) + '\t' + score);
						} else if (variant.isPhased(founderMember) && variant.isPhased(QuartetMember.KID1) && variant.isPhased(QuartetMember.KID2) && variant.isHeterozygous(founderMember)) {						
							int score;
							if (variant.getAlleles(QuartetMember.KID1)[founderAllele] == variant.getAlleles(QuartetMember.KID2)[founderAllele]) {
								score = 1;
							} else {
								score = -1;
							}
							System.out.println(variant.getChromosome() + '\t' + variant.getPosition() + '\t' + (variant.getPosition() + 1) + '\t' + score);
						}
					} catch (InvalidVCFLineException | FilteredVCFLineException | InvalidVCFFieldException | PartiallyCalledVariantException e) {
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
