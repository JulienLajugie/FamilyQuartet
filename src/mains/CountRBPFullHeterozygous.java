package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import dataStructures.PhasedGenotypesSeries;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;
import dataStructures.QuartetMember;
import dataStructures.SegmentalDuplication;
import dataStructures.SegmentalDuplicationList;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Counts the number of full heterozygous variants phased using the RBP algorithm 
 * @author Julien Lajugie
 */
public class CountRBPFullHeterozygous {

	/**
	 * Usage: java CountRBPFullHeterozygous.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>
	 * @param args -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) {
			System.out.println("Usage: java CountRBPFullHeterozygous.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>");
			System.exit(-1);
		} else {
			File geneticPhasingFile = null;
			File physicalPhasingFile = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-g")) {
					geneticPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					physicalPhasingFile = new File(args[i + 1]);
				}
			}
			try {
				countRBPFullHeterozygous(geneticPhasingFile, physicalPhasingFile);
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
		// case with no -g parameter
		if (!args[0].equals("-g") && !args[2].equals("-g")) {
			return false;
		}
		// case with no -p parameter
		if (!args[0].equals("-p") && !args[2].equals("-p")) {
			return false;
		}
		return true;
	}


	/**
	 * Counts the number of full heterozygous variants phased using the RBP algorithm
	 * @param geneticPhasingFile file containing the result of the genetic phasing
	 * @param physicalPhasingFile file containing the result of the physical phasing
	 * @throws IOException
	 */
	private static void countRBPFullHeterozygous(File geneticPhasingFile, File physicalPhasingFile) throws IOException {
		BufferedReader reader = null;
		String line = null;
		int fullHeterozygousPhasedCount = 0;
		int fullHeterozygousUnphasedCount = 0;
		Map<QuartetMember, SegmentalDuplicationList> commonPhasedBlocks = createCommonPhasedBlocks(geneticPhasingFile, physicalPhasingFile);
		Map<QuartetMember, SegmentalDuplicationList> RBPhasedBlocks = createRBPBlocks(physicalPhasingFile);
		try {
			reader = new BufferedReader(new FileReader(physicalPhasingFile));
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant variant;
						variant = new Variant(line);
						if (variant.getGenotypePattern().equals("ab/ab;ab/ab")) {
							boolean phasedMemberFound = false;
							for (QuartetMember member: QuartetMember.values()) {
								SegmentalDuplication RBPBlock = RBPhasedBlocks.get(member).getBlock(variant.getChromosome(), variant.getPosition());
								if (RBPBlock != null) {
									if (commonPhasedBlocks.get(member).containsBlockOverlapping(variant.getChromosome(), RBPBlock)) {
										phasedMemberFound = true;
									}
								}
							}
							// add one to the vector count if a member is phased (the other members can be deducted using the inheritance state blocks)
							if (phasedMemberFound) {
								fullHeterozygousPhasedCount++;
							} else {
								fullHeterozygousUnphasedCount++;								
							}
						}
					} catch (VCFException e) {} 					
				}				
			}
			System.out.println("Full heterozygous vectors phased: " + fullHeterozygousPhasedCount);
			System.out.println("Full heterozygous vectors unphased: " + fullHeterozygousUnphasedCount);
		} finally {
			if (reader != null) {
				reader.close();
			}
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
