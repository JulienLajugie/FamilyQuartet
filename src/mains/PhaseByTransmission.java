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
import dataStructures.QuartetMember;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Phases a quartet using a transmission algorithm
 * @author Julien Lajugie
 */
public class PhaseByTransmission {

	/**
	 * Usage: java PhaseByTransmission.java -b <path to the block bgr file> -v <path to the transmission phased vcf file>
	 * @param args -b <path to the block bgr file> -v <path to the transmission phased vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java PhaseByTransmission.java -b <path to the block bgr file> -v <path to the transmission phased vcf file>");
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
	 * Phases a quartet using a transmission algorithm
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
						// full heterozygous vectors can't be phased
						if (!variant.getGenotypePattern().equals("ab/ab;ab/ab")) {
							// first we phase every homozygous variant
							for (QuartetMember member: QuartetMember.values()) {
								if (variant.isHomozygous(member)) {
									variant.setPhase(member, true);
								}
							}
							InheritanceStateBlock<CrossTriosInheritanceState> block = blockList.getBlock(variant);
							// cannot phase if there is no block information
							CrossTriosInheritanceState blockState = null;
							if (block != null) {
								blockState = block.getBlockState();
							}
							// cannot phase MIEs and SCEs
							if (!variant.isMIE() && ((blockState == null) || !variant.isSCE(blockState))) {
								if (variant.isHomozygous(QuartetMember.FATHER)) {
									phaseFromHomozygousFather(variant);
								} else if (variant.isHomozygous(QuartetMember.MOTHER)) {
									phaseFromHomozygousMother(variant);
								} else if (variant.isHomozygous(QuartetMember.KID1)) {
									phaseFromHomozygousKid1(variant, blockState);
								} else if (variant.isHomozygous(QuartetMember.KID2)) {
									phaseFromHomozygousKid2(variant, blockState);
								}
							}
						}
						lineToPrint = substituteVcfLine(line, variant);
					} catch (VCFException exception) {}
				}
				System.out.println(lineToPrint);
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Phase the quartet using the father genotype
	 * @param variant variant to phase
	 */
	private static void phaseFromHomozygousFather(Variant variant) {
		AlleleType fatherAlleles = variant.getAlleles(QuartetMember.FATHER)[0];
		AlleleType otherAllele = fatherAlleles.getOpositeAllele();
		if (!variant.isHomozygous(QuartetMember.KID1)) {
			variant.setGenotype(QuartetMember.KID1, otherAllele, fatherAlleles, true);
		}
		if (!variant.isHomozygous(QuartetMember.KID2)) {
			variant.setGenotype(QuartetMember.KID2, otherAllele, fatherAlleles, true);
		}
		if (!variant.isHomozygous(QuartetMember.MOTHER)) {
			AlleleType kid1MaternalAllele = variant.getAlleles(QuartetMember.KID1)[0];
			variant.setGenotype(QuartetMember.MOTHER, kid1MaternalAllele, kid1MaternalAllele.getOpositeAllele(), true);
		}
	}


	/**
	 * Phase the quartet using the mother genotype
	 * @param variant variant to phase
	 * @param blockState state of the block to phase
	 */
	private static void phaseFromHomozygousMother(Variant variant) {
		AlleleType motherAlleles = variant.getAlleles(QuartetMember.MOTHER)[0];
		AlleleType otherAllele = motherAlleles.getOpositeAllele();
		if (!variant.isHomozygous(QuartetMember.KID1)) {
			variant.setGenotype(QuartetMember.KID1, motherAlleles, otherAllele, true);
		}
		if (!variant.isHomozygous(QuartetMember.KID2)) {
			variant.setGenotype(QuartetMember.KID2, motherAlleles, otherAllele, true);
		}
		if (!variant.isHomozygous(QuartetMember.FATHER)) {
			AlleleType kid1PaternalAllele = variant.getAlleles(QuartetMember.KID1)[1];
			variant.setGenotype(QuartetMember.FATHER, kid1PaternalAllele, kid1PaternalAllele.getOpositeAllele(), true);
		}
	}	

	/**
	 * Phase the quartet using the kid1 genotype
	 * @param variant variant to phase
	 * @param blockState state of the block to phase
	 */
	private static void phaseFromHomozygousKid1(Variant variant, CrossTriosInheritanceState blockState) {
		AlleleType kid1Alleles = variant.getAlleles(QuartetMember.KID1)[0];
		AlleleType otherAllele = kid1Alleles.getOpositeAllele();
		if (!variant.isHomozygous(QuartetMember.FATHER)) {
			variant.setGenotype(QuartetMember.FATHER, kid1Alleles, otherAllele, true);
		}
		if (!variant.isHomozygous(QuartetMember.MOTHER)) {
			variant.setGenotype(QuartetMember.MOTHER, kid1Alleles, otherAllele, true);
		}
		if (blockState != null) {
			if (!variant.isHomozygous(QuartetMember.KID2)) {
				if (blockState.getPaternalTrioState() != TrioInheritanceState.UNKNOWN) {
					if (blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) {
						variant.setGenotype(QuartetMember.KID2, otherAllele, kid1Alleles, true);
					} else if (blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
						variant.setGenotype(QuartetMember.KID2, kid1Alleles, otherAllele, true);
					}
				} else if (blockState.getMaternalTrioState() != TrioInheritanceState.UNKNOWN) {
					if (blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) {
						variant.setGenotype(QuartetMember.KID2, kid1Alleles, otherAllele, true);
					} else if (blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
						variant.setGenotype(QuartetMember.KID2, otherAllele, kid1Alleles, true);
					}
				}
			}
		}
	}


	/**
	 * Phase the quartet using the kid2 genotype
	 * @param variant variant to phase
	 * @param blockState state of the block to phase
	 */
	private static void phaseFromHomozygousKid2(Variant variant, CrossTriosInheritanceState blockState) {
		if (blockState != null) {
			AlleleType kid2Alleles = variant.getAlleles(QuartetMember.KID2)[0];
			AlleleType otherAllele = kid2Alleles.getOpositeAllele();
			if (blockState.getPaternalTrioState() != TrioInheritanceState.UNKNOWN) {
				if (blockState.getPaternalTrioState() == TrioInheritanceState.IDENTICAL) {
					// case paternal identical
					if (!variant.isHomozygous(QuartetMember.FATHER)) {
						variant.setGenotype(QuartetMember.FATHER, kid2Alleles, otherAllele, true);
					}
					if (!variant.isHomozygous(QuartetMember.MOTHER)) {
						variant.setGenotype(QuartetMember.MOTHER, otherAllele, kid2Alleles, true);
					}
				} else if (blockState.getPaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
					// case paternal non-identical
					if (!variant.isHomozygous(QuartetMember.FATHER)) {
						variant.setGenotype(QuartetMember.FATHER, otherAllele, kid2Alleles, true);
					}
					if (!variant.isHomozygous(QuartetMember.MOTHER)) {
						variant.setGenotype(QuartetMember.MOTHER, kid2Alleles, otherAllele, true);
					}
				}
			} else if (blockState.getMaternalTrioState() != TrioInheritanceState.UNKNOWN) {
				if (blockState.getMaternalTrioState() == TrioInheritanceState.IDENTICAL) {
					// case maternal identical
					if (!variant.isHomozygous(QuartetMember.FATHER)) {
						variant.setGenotype(QuartetMember.FATHER, otherAllele, kid2Alleles, true);
					}
					if (!variant.isHomozygous(QuartetMember.MOTHER)) {
						variant.setGenotype(QuartetMember.MOTHER, kid2Alleles, otherAllele, true);
					}
				} else if (blockState.getMaternalTrioState() == TrioInheritanceState.NON_IDENTICAL) {
					// case paternal non-identical
					if (!variant.isHomozygous(QuartetMember.FATHER)) {
						variant.setGenotype(QuartetMember.FATHER, kid2Alleles, otherAllele, true);
					}
					if (!variant.isHomozygous(QuartetMember.MOTHER)) {
						variant.setGenotype(QuartetMember.MOTHER, otherAllele, kid2Alleles, true);
					}
				}
			}
			AlleleType maternalKid1Allele = variant.getAlleles(QuartetMember.MOTHER)[0];
			AlleleType paternalKid1Allele = variant.getAlleles(QuartetMember.FATHER)[0];
			variant.setGenotype(QuartetMember.KID1, maternalKid1Allele, paternalKid1Allele, true);
		}
	}


	/**
	 * @param line input line
	 * @param variant phased variant
	 * @return the input line with the unphased genotype substituted by a phase genotype 
	 */
	private static String substituteVcfLine(String line, Variant variant) {
		String[] splitLine = line.split("\t");
		String phasedVcfLine = "";
		int i;
		for (i = 0; i < 9; i++) {
			phasedVcfLine += splitLine[i] + "\t";
		}
		for (QuartetMember member: QuartetMember.values()) {
			String memberPhasing = variant.isPhased(member) == true ? "|" : "/";
			String phasedMemberGenotype = variant.getAlleles(member)[0].getIntValue() + memberPhasing + variant.getAlleles(member)[1].getIntValue();
			String memberGenotypeField = phasedMemberGenotype + splitLine[i++].trim().substring(3); 
			phasedVcfLine += memberGenotypeField + "\t";
		}
		return phasedVcfLine;
	}
}
