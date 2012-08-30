package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.AlleleType;
import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.QuartetInheritanceState;
import dataStructures.QuartetMember;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Improves the result of the double transmission phasing by phasing the variants 
 * that are heterozygous in both parents as well as in one and only one of the kids.
 * These variant can't be phased in the context of a trio but can't be phased using the 
 * inheritance genotype block information.
 * @author Julien Lajugie
 */
public class ExtendDoubleTransmissionPhasing {

	/**
	 * Usage: java ExtendDoubleTransmissionPhasing.java -b <path to the block bgr file> -v <path to the transmission phased vcf file>
	 * @param args -b <path to the block bgr file> -v <path to the transmission phased vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ExtendDoubleTransmissionPhasing.java -b <path to the block bgr file> -v <path to the transmission phased vcf file>");
			System.exit(-1);
		} else {
			File blockFile = null;
			File vcfFile = null;
			for (int i = 0; i <= 2; i += 2) {
				if (args[i].equals("-b")) {
					blockFile = new File(args[i + 1]);
				}
				if (args[i].equals("-v")) {
					vcfFile = new File(args[i + 1]);
				}
			}
			try {
				extendDoubleTransmissionPhasing(blockFile, vcfFile);
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
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b")) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v")) {
			return false;
		}
		return true;
	}



	/**
	 * Improves the result of the double transmission phasing by phasing the variants 
	 * that are heterozygous in both parents as well as in one and only one of the kids.
	 * These variant can't be phased in the context of a trio but can't be phased using the 
	 * inheritance genotype block information
	 * @param blockFile bedgraph with the blocks
	 * @param vcfFile vcf file with the result of the two trios phasing using the transmission method
	 * @throws IOException
	 */
	private static void extendDoubleTransmissionPhasing(File blockFile, File vcfFile) throws IOException {
		// load the block file
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList;
		blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		BufferedReader reader = null;
		 //int phasableVariantCount = 0;
		 //int phasableSCEVariantCount = 0;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				String lineToPrint = line;
				if (line.charAt(0) != '#') {
					try {
						Variant variant = new Variant(line);
						// we work only with the variants heterozygote for both parents and one child (ie with a "ab/ab;aa/ab" genotype pattern)
						if (variant.getGenotypePattern().equals("ab/ab;aa/ab")) {
							InheritanceStateBlock<CrossTriosInheritanceState> block = blockList.getBlock(variant);
							if (block != null) {
								//phasableVariantCount++;
								CrossTriosInheritanceState blockState = blockList.getBlock(variant).getBlockState();
								if (!blockState.isCompatibleWith(QuartetInheritanceState.MATERNAL) && !blockState.isCompatibleWith(QuartetInheritanceState.PATERNAL)) {
									//phasableSCEVariantCount++;
								} else {
									if (variant.isHomozygous(QuartetMember.KID1)) {								
										phaseKid2(variant, blockState);
									} else if (variant.isHomozygous(QuartetMember.KID2)){	
										phaseTrioParentsKid1(variant, blockState);
									}
									lineToPrint = substituteVcfLine(line, variant);
								}
							}
						}
					} catch (VCFException exception) {
						// do nothing
					}
				}
				System.out.println(lineToPrint);
			}
			 //System.out.println("Number of variants with \"ab/ab;aa/ab\" genotype pattern: " + phasableVariantCount);
			 //System.out.println("% of SCE within this subset: " + (phasableSCEVariantCount / (double) phasableVariantCount));
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Phase the 1st child and the parents of the quartet ("ab/ab;ab+aa")
	 * @param variant variant to phase
	 * @param blockState state of the block to phase
	 */
	private static void phaseTrioParentsKid1(Variant variant, CrossTriosInheritanceState blockState) {		
		AlleleType child2Alleles = variant.getAlleles(QuartetMember.KID2)[0];
		AlleleType otherAllele = child2Alleles.getOpositeAllele();
		if (blockState.getPaternalTrioState() != TrioInheritanceState.UNKNOWN) {
			if (blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) {
				// case paternal identical
				variant.setGenotype(QuartetMember.FATHER, child2Alleles, otherAllele, true);
				variant.setGenotype(QuartetMember.MOTHER, otherAllele, child2Alleles, true);
				variant.setGenotype(QuartetMember.KID1, otherAllele, child2Alleles, true);				
			} else if (blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
				// case paternal non-identical
				variant.setGenotype(QuartetMember.FATHER, otherAllele, child2Alleles, true);
				variant.setGenotype(QuartetMember.MOTHER, child2Alleles, otherAllele, true);
				variant.setGenotype(QuartetMember.KID1, child2Alleles, otherAllele, true);
			}
		} else if (blockState.getMaternalTrioState() != TrioInheritanceState.UNKNOWN) {
			if (blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) {
				// case maternal identical
				variant.setGenotype(QuartetMember.FATHER, otherAllele, child2Alleles, true);
				variant.setGenotype(QuartetMember.MOTHER, child2Alleles, otherAllele, true);
				variant.setGenotype(QuartetMember.KID1, child2Alleles, otherAllele, true);
			} else if (blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
				// case paternal non-identical
				variant.setGenotype(QuartetMember.FATHER, child2Alleles, otherAllele, true);
				variant.setGenotype(QuartetMember.MOTHER, otherAllele, child2Alleles, true);
				variant.setGenotype(QuartetMember.KID1, otherAllele, child2Alleles, true);			
			}
		}
	}


	/**
	 * Phase the second kid of the quartet ("ab/ab;aa+ab" case)
	 * @param variant variant to phase
	 * @param blockState state of the block
	 */
	private static void phaseKid2(Variant variant, CrossTriosInheritanceState blockState) {
		if (blockState.getPaternalTrioState() != TrioInheritanceState.UNKNOWN) {
			if (blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) {
				// case paternal identical
				variant.setGenotype(QuartetMember.KID2, variant.getAlleles(QuartetMember.FATHER)[1], variant.getAlleles(QuartetMember.FATHER)[0], true);
			} else if (blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
				// case paternal non-identical
				variant.setGenotype(QuartetMember.KID2, variant.getAlleles(QuartetMember.FATHER)[0], variant.getAlleles(QuartetMember.FATHER)[1], true);
			}
		} else if (blockState.getMaternalTrioState() != TrioInheritanceState.UNKNOWN) {
			// if the father state is unknown we use the mother state for the phasing
			if (blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) {
				// case maternal identical
				variant.setGenotype(QuartetMember.KID2, variant.getAlleles(QuartetMember.MOTHER)[0], variant.getAlleles(QuartetMember.MOTHER)[1], true);
			} else if (blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
				// case maternal non-identical
				variant.setGenotype(QuartetMember.KID2, variant.getAlleles(QuartetMember.MOTHER)[1], variant.getAlleles(QuartetMember.MOTHER)[0], true);
			}			
		}
	}
	
	
	/**
	 * @param line input line
	 * @param variant phased variant
	 * @return the input line with the unphased genotype substituted by a phase genotype 
	 */
	private static String substituteVcfLine(String line, Variant variant) {
		String[] splitLine = line.split("\t");
		
		String phasedFatherGenotype = variant.getAlleles(QuartetMember.FATHER)[0].getIntValue() + "|" + variant.getAlleles(QuartetMember.FATHER)[1].getIntValue(); 
		String phasedmotherGenotype = variant.getAlleles(QuartetMember.MOTHER)[0].getIntValue() + "|" + variant.getAlleles(QuartetMember.MOTHER)[1].getIntValue(); 
		String phasedkid1Genotype = variant.getAlleles(QuartetMember.KID1)[0].getIntValue() + "|" + variant.getAlleles(QuartetMember.KID1)[1].getIntValue(); 
		String phasedkid2Genotype = variant.getAlleles(QuartetMember.KID2)[0].getIntValue() + "|" + variant.getAlleles(QuartetMember.KID2)[1].getIntValue();
		
		String fatherGenotype = phasedFatherGenotype + splitLine[9].trim().substring(3);
		String motherGenotype = phasedmotherGenotype + splitLine[10].trim().substring(3);
		String kid1Genotype = phasedkid1Genotype + splitLine[11].trim().substring(3);
		String kid2Genotype = phasedkid2Genotype + splitLine[12].trim().substring(3);
		
		String phasedVcfLine = "";
		for (int i = 0; i < 9; i++) {
			phasedVcfLine += splitLine[i] + "\t";
		}
		phasedVcfLine += fatherGenotype + "\t";
		phasedVcfLine += motherGenotype + "\t";
		phasedVcfLine += kid1Genotype + "\t";
		phasedVcfLine += kid2Genotype;
		return phasedVcfLine;
	}
}
