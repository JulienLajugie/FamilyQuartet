package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.InheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;

public class GenerateVCFForHaploscript {


	/**
	 * Usage: java GenerateVCFForHaploscript -v <path to the VCF file> -b <path to the block file> -c <path to the compression block>
	 * @param args -v <path to the VCF file> -b <path to the block file> -c <path to the compression block>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateVCFForHaploscript -v <path to the VCF file> -b <path to the block file> -c <path to the compression block>");
			System.exit(-1);
		} else {
			File VCFFile = null, blockFile = null, compressionFile = null;
			for (int i = 0; i <= 4; i += 2) {
				if (args[i].equals("-v")) {
					VCFFile = new File(args[i + 1]);
				}
				if (args[i].equals("-b")) {
					blockFile = new File(args[i + 1]);
				}
				if (args[i].equals("-c")) {
					compressionFile = new File(args[i + 1]);
				}
			}
			try {
				generateVCFForHaploscript(VCFFile, blockFile, compressionFile);
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
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v") && !args[4].equals("-v")) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b")) {
			return false;
		}
		// case with no -c parameter
		if (!args[0].equals("-c") && !args[2].equals("-c") && !args[4].equals("-c")) {
			return false;
		}
		return true;

	}


	/**
	 * Generates a VCF file filtered for Haplotyping
	 * @param VCFFile input VCF file
	 * @param blockFile ISCA file with the inheritance state blocks
	 * @param compressionFile bgr file with the compression blocks 
	 * @throws IOException
	 */
	private static void generateVCFForHaploscript(File VCFFile, File blockFile, File compressionFile) throws IOException {
		// load the inheritance block list
		InheritanceStateBlockList blockList = new InheritanceStateBlockList();
		blockList.loadFromISCAFile(blockFile);
		// load the compression block list
		InheritanceStateBlockList compressionList = new InheritanceStateBlockList();
		blockList.loadFromBgrFile(compressionFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			String previousLine = null;
			int previousVariantPosition = -1;
			String previousChromo = "";
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process variants with more than one alternative allele or indels
						if ((currentVariant.getAlternatievAllele().length() != 1) || (currentVariant.getReferenceAllele().length() != 1)) {
							throw new InvalidVCFLineException();
						}
						// we don't want MIE
						if ((currentVariant.getInheritanceStates()[0] == InheritanceState.MIE)) {
							throw new InvalidVCFLineException();
						}
						// we don't want SCE
						InheritanceStateBlock variantInheritanceBlock = blockList.getBlock(currentVariant);
						if ((variantInheritanceBlock != null) && (variantInheritanceBlock.isSCE(currentVariant))) {
							throw new InvalidVCFLineException();
						}
						// we don't want variants in compression blocks
						InheritanceStateBlock variantCompressionBlock = compressionList.getBlock(currentVariant);
						if (variantCompressionBlock != null) {
							throw new InvalidVCFLineException();
						}						
						// we don't want to have two vcf lines for the same position
						if ((currentVariant.getPosition() == previousVariantPosition) &&
								(currentVariant.getChromosome().equals(previousChromo))) {
							previousLine = null;
						} else {							
							if (previousLine != null) {
								System.out.println(previousLine);
							}
							previousChromo = currentVariant.getChromosome();
							previousVariantPosition = currentVariant.getPosition();
							previousLine = line;
						}
						
					} catch (InvalidVCFLineException e) {
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
