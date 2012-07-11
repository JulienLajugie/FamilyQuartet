package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.CrossTriosInheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.PhasedGenotypesSeries;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;
import dataStructures.QuartetMember;
import dataStructures.SegmentalDuplication;
import dataStructures.SegmentalDuplicationList;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Phases the full heterozygous variants of the TP using the result of the RBP 
 * @author Julien Lajugie
 */
public class PhaseFullHeterozygous {

	/**
	 * Usage: java PhaseFullHeterozygous.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -b <path to the inheritance block file>
	 * @param args -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -b <path to the inheritance block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) {
			System.out.println("Usage: java PhaseFullHeterozygous.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -b <path to the inheritance block file>");
			System.exit(-1);
		} else {
			File geneticPhasingFile = null;
			File physicalPhasingFile = null;
			File inheritanceBlockFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-g")) {
					geneticPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					physicalPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-b")) {
					inheritanceBlockFile = new File(args[i + 1]);
				}				
			}
			try {
				countRBPFullHeterozygous(geneticPhasingFile, physicalPhasingFile, inheritanceBlockFile);
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
		// case with no -g parameter
		if (!args[0].equals("-g") && !args[2].equals("-g") && !args[4].equals("-g")) {
			return false;
		}
		// case with no -p parameter
		if (!args[0].equals("-p") && !args[2].equals("-p") && !args[4].equals("-p")) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b")) {
			return false;
		}
		return true;
	}


	/**
	 * Counts the number of full heterozygous variants phased using the RBP algorithm
	 * @param geneticPhasingFile file containing the result of the genetic phasing
	 * @param physicalPhasingFile file containing the result of the physical phasing
	 * @param inheritanceBlockFile file containing the inheritance state blocks
	 * @throws IOException
	 */
	private static void countRBPFullHeterozygous(File geneticPhasingFile, File physicalPhasingFile, File inheritanceBlockFile) throws IOException {
		BufferedReader reader = null;
		String line = null;
		Map<QuartetMember, SegmentalDuplicationList> commonPhasedBlocks = createCommonPhasedBlocks(geneticPhasingFile, physicalPhasingFile);
		Map<QuartetMember, SegmentalDuplicationList> RBPhasedBlocks = createRBPBlocks(physicalPhasingFile);
		PhasedVectorList RBPVectors = new PhasedVectorList();
		RBPVectors.loadFromVCFFile(physicalPhasingFile);
		PhasedVectorList TPVectors = new PhasedVectorList();
		TPVectors.loadFromVCFFile(geneticPhasingFile);
		InheritanceStateBlockList<CrossTriosInheritanceState> isBlockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(inheritanceBlockFile);		

		try {
			reader = new BufferedReader(new FileReader(geneticPhasingFile));
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					System.out.println(line);
				} else {
					try {
						Variant variant = new Variant(line);
						String chromosome = variant.getChromosome();
						CrossTriosInheritanceStateBlock isBlock = (CrossTriosInheritanceStateBlock) isBlockList.getBlock(variant); 
						// indels are not phased by the RBP software
						if (!variant.isIndel() && !variant.isMIE() && ((isBlock.getBlockState() == null) || !variant.isSCE(isBlock.getBlockState()))) {
							if (variant.getGenotypePattern().equals("ab/ab;ab/ab")) {
								boolean isVariantPhased = false;
								PhasedVector variantVector = RBPVectors.getPhasedVector(chromosome, variant.getPosition());
								// we unphase the RBP vector
								for (QuartetMember member: QuartetMember.values()) {
									variantVector.setPhasing(member, false);
								}
								// we try to phase it
								for (QuartetMember member: QuartetMember.values()) {
									if (!isVariantPhased) {
										SegmentalDuplication RBPBlock = RBPhasedBlocks.get(member).getBlock(variant.getChromosome(), variant.getPosition());
										if (RBPBlock != null) {
											SegmentalDuplication commonBlock = commonPhasedBlocks.get(member).getBlockOverlapping(variant.getChromosome(), RBPBlock); 
											int commonVariantPosition = commonBlock.getStartPosition();
											if (commonBlock != null) {
												String RBPGenotype = RBPVectors.getPhasedVector(chromosome, commonVariantPosition).getGenotype(member);
												String TPGenotype = TPVectors.getPhasedVector(chromosome, commonVariantPosition).getGenotype(member);
												boolean needToBeInverted = RBPGenotype.equals(TPGenotype);
												if (needToBeInverted) {
													variantVector.invert(member);
												}
												variantVector.setPhasing(member, true);
												if (isBlock != null) {
													isVariantPhased = true;
													phaseFamily(member, variantVector, isBlock.getBlockState());
												}
											}
										}
									}
								}
								line = phaseVCFLine(line, variantVector);
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
	 * Phases the rest of the familly using one phased member and the inheritance state information
	 * @param phasedMember member that is phased
	 * @param vectorToPhase the vector to phase
	 * @param inheritanceState the inheritance state of the vector
	 */
	private static void phaseFamily(QuartetMember phasedMember, PhasedVector vectorToPhase, CrossTriosInheritanceState inheritanceState) {
		TrioInheritanceState maternalState = inheritanceState.getMaternalTrioState();
		if (maternalState == TrioInheritanceState.UNKNOWN) {
			// a ab/ab;ab/ab variant not SCE can only be identical or non-identical which implies that the maternal and paternal states are indentical
			maternalState = inheritanceState.getPaternalTrioState();
		}
		assert (maternalState == TrioInheritanceState.IDENTICAL) || (maternalState == TrioInheritanceState.NON_IDENTICAL);
		// we start by phasing the kid1
		switch (phasedMember) {
		case FATHER:
			if (vectorToPhase.getGenotype(QuartetMember.FATHER).charAt(0) != vectorToPhase.getGenotype(QuartetMember.KID1).charAt(1)) {
				vectorToPhase.invert(QuartetMember.KID1);
			}
			vectorToPhase.setPhasing(QuartetMember.KID1, true);
			break;
		case MOTHER:
			if (vectorToPhase.getGenotype(QuartetMember.MOTHER).charAt(0) != vectorToPhase.getGenotype(QuartetMember.KID1).charAt(0)) {
				vectorToPhase.invert(QuartetMember.KID1);
			}
			vectorToPhase.setPhasing(QuartetMember.KID1, true);
			break;
		case KID1:
			// nothing to do
			break;
		case KID2:
			char kid1Allele1 = vectorToPhase.getGenotype(QuartetMember.KID1).charAt(0);
			char kid2Allele1 = vectorToPhase.getGenotype(QuartetMember.KID2).charAt(0);
			if (((maternalState == TrioInheritanceState.IDENTICAL) && (kid1Allele1 != kid2Allele1)) ||
					((maternalState == TrioInheritanceState.NON_IDENTICAL) && (kid1Allele1 == kid2Allele1))) {
				vectorToPhase.invert(QuartetMember.KID1);				
			}
			vectorToPhase.setPhasing(QuartetMember.KID1, true);
			break;			
		}
		// we now phase the remaining of the familly
		if (!vectorToPhase.isPhased(QuartetMember.FATHER)) {
			if (vectorToPhase.getGenotype(QuartetMember.FATHER).charAt(0) != vectorToPhase.getGenotype(QuartetMember.KID1).charAt(1)) {
				vectorToPhase.invert(QuartetMember.FATHER);
			}
			vectorToPhase.setPhasing(QuartetMember.FATHER, true);
		}
		if (!vectorToPhase.isPhased(QuartetMember.MOTHER)) {
			if (vectorToPhase.getGenotype(QuartetMember.MOTHER).charAt(0) != vectorToPhase.getGenotype(QuartetMember.KID1).charAt(0)) {
				vectorToPhase.invert(QuartetMember.MOTHER);
			}
			vectorToPhase.setPhasing(QuartetMember.MOTHER, true);
		}
		if (!vectorToPhase.isPhased(QuartetMember.KID2)) {
			char kid1Allele1 = vectorToPhase.getGenotype(QuartetMember.KID1).charAt(0);
			char kid2Allele1 = vectorToPhase.getGenotype(QuartetMember.KID2).charAt(0);
			if (((maternalState == TrioInheritanceState.IDENTICAL) && (kid1Allele1 != kid2Allele1)) ||
					((maternalState == TrioInheritanceState.NON_IDENTICAL) && (kid1Allele1 == kid2Allele1))) {
				vectorToPhase.invert(QuartetMember.KID2);				
			}
			vectorToPhase.setPhasing(QuartetMember.KID2, true);
		}
	}


	/**
	 * @param line unphased VCF line
	 * @param phasedVector vector with phased genotypes
	 * @return a vcf line where the unphased genotypes that can be phased using the vector are phased
	 */
	private static String phaseVCFLine(String line, PhasedVector phasedVector) {
		String[] splitLine = line.split("\t");
		for (QuartetMember member: QuartetMember.values()) {
			if (phasedVector.isPhased(member)) {
				String phasedGenotype = phasedVector.getGenotype(member);
				int memberInfoFieldIndex = getMemberInfoFieldIndex(member);
				String infoFieldWithoutGenotype = splitLine[memberInfoFieldIndex].trim().substring(3);
				splitLine[memberInfoFieldIndex] = phasedGenotype + infoFieldWithoutGenotype;
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


	/**
	 * @param physicalPhasingFile vcf file phased using a physical algorithm
	 * @return a map with the phased blocks starting and ending by an heterozygous variant for each family member
	 * @throws IOException
	 */
	private static Map<QuartetMember, SegmentalDuplicationList> createRBPBlocks(File physicalPhasingFile) throws IOException {
		// create map with genetic and RBP phased block lists
		Map<QuartetMember, SegmentalDuplicationList> phasedBlockLists = new HashMap<>();
		phasedBlockLists.put(QuartetMember.FATHER, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.MOTHER, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.KID1, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.KID2, new SegmentalDuplicationList());

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(physicalPhasingFile));
			String line = null;

			// start positions
			Map<QuartetMember, Integer> startPositions = new HashMap<>();
			startPositions.put(QuartetMember.FATHER, 0);
			startPositions.put(QuartetMember.MOTHER, 0);
			startPositions.put(QuartetMember.KID1, 0);
			startPositions.put(QuartetMember.KID2, 0);

			// stop positions
			Map<QuartetMember, Integer> stopPositions = new HashMap<>();
			stopPositions.put(QuartetMember.FATHER, 0);
			stopPositions.put(QuartetMember.MOTHER, 0);
			stopPositions.put(QuartetMember.KID1, 0);
			stopPositions.put(QuartetMember.KID2, 0);

			// chromosomes
			Map<QuartetMember, String> chromosomes = new HashMap<>();
			chromosomes.put(QuartetMember.FATHER, null);
			chromosomes.put(QuartetMember.MOTHER, null);
			chromosomes.put(QuartetMember.KID1, null);
			chromosomes.put(QuartetMember.KID2, null);

			// true if the variant for the specified member is the first of the block
			Map<QuartetMember, Boolean> isBlockFirstVariants = new HashMap<>();
			isBlockFirstVariants.put(QuartetMember.FATHER, true);
			isBlockFirstVariants.put(QuartetMember.MOTHER, true);
			isBlockFirstVariants.put(QuartetMember.KID1, true);
			isBlockFirstVariants.put(QuartetMember.KID2, true);

			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					Variant currentVariant;
					try {
						currentVariant = new Variant(line);
						for (QuartetMember member: QuartetMember.values()) {
							// case where the variant is not phased or is on a new chromosome (meaning that the previous block ended)
							if ((!currentVariant.isPhased(member)) || (!currentVariant.getChromosome().equals(chromosomes.get(member)))) {
								// case where the previous block is not empty
								if (!isBlockFirstVariants.get(member)) {
									SegmentalDuplication dupToAdd = new SegmentalDuplication(startPositions.get(member), stopPositions.get(member));
									phasedBlockLists.get(member).addDuplication(chromosomes.get(member), dupToAdd);
									isBlockFirstVariants.put(member, true);
								}
								chromosomes.put(member, currentVariant.getChromosome());
								startPositions.put(member, currentVariant.getPosition());
							}
							// case where the variant is phased with the previsous one
							if ((currentVariant.isPhased(member)) && (currentVariant.isHeterozygous(member))) {
								isBlockFirstVariants.put(member, false);
								stopPositions.put(member, currentVariant.getPosition());
							} 
						}
					} catch (VCFException e) {
						// do nothing
					}
				}
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
		return phasedBlockLists;
	}


	/**
	 * 
	 * @param geneticPhasingFile vcf file phased using a genetic algorithm
	 * @param physicalPhasingFile vcf file phased using a physical algorithm
	 * @return a map with the phased blocks (with only vector phased for both genetic and physical method) variant for each family member
	 * @throws IOException
	 */
	private static Map<QuartetMember, SegmentalDuplicationList> createCommonPhasedBlocks(File geneticPhasingFile, File physicalPhasingFile) throws IOException {
		// load genetic phasing file		
		PhasedVectorList geneticVectorList = new PhasedVectorList();
		geneticVectorList.loadFromVCFFile(geneticPhasingFile);

		// load physical phasing file
		PhasedVectorList physicalVectorList = new PhasedVectorList();
		physicalVectorList.loadFromVCFFile(physicalPhasingFile);

		// create map with the phased series
		Map<QuartetMember, PhasedGenotypesSeries> phasedSeries = new HashMap<>();
		phasedSeries.put(QuartetMember.FATHER, new PhasedGenotypesSeries());
		phasedSeries.put(QuartetMember.MOTHER, new PhasedGenotypesSeries());
		phasedSeries.put(QuartetMember.KID1, new PhasedGenotypesSeries());
		phasedSeries.put(QuartetMember.KID2, new PhasedGenotypesSeries());

		// create map with genetic and RBP phased block lists
		Map<QuartetMember, SegmentalDuplicationList> phasedBlockLists = new HashMap<>();
		phasedBlockLists.put(QuartetMember.FATHER, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.MOTHER, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.KID1, new SegmentalDuplicationList());
		phasedBlockLists.put(QuartetMember.KID2, new SegmentalDuplicationList());

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(physicalPhasingFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					PhasedVector geneticPhasedVector = geneticVectorList.getPhasedVector(chromosome, position);
					PhasedVector physicalPhasedVector = physicalVectorList.getPhasedVector(chromosome, position);
					if ((geneticPhasedVector != null) && (physicalPhasedVector != null)) {
						for (QuartetMember member: QuartetMember.values()) {
							int result = phasedSeries.get(member).addGeneticPhysicalGenotypes(geneticPhasedVector.getGenotype(member), physicalPhasedVector.getGenotype(member), chromosome, position);
							if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
								SegmentalDuplication phasedBlock = phasedSeries.get(member).getBlock();
								if (phasedBlock != null) {
									phasedBlockLists.get(member).addDuplication(chromosome, phasedBlock);
								}
								phasedSeries.get(member).reset();
							}
						}
					}
				}
			}
			for (QuartetMember member: QuartetMember.values()) {
				phasedBlockLists.get(member).sort();
			}
			return phasedBlockLists;
		} finally {
			if (reader != null) {
				reader.close();
			}
		}	
	}
}
