package mains;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;
import dataStructures.PhasedVectorsInsert;
import dataStructures.QuartetMember;

/**
 * Tries to improve the result of a genetic phasing (eg: haplotyping, transmission phasing) 
 * using the result of a physical phasing (ie read back phasing). 
 * @author Julien Lajugie
 */
public class ExtendPhasingUsingRBP {


	private final static String[] CHROMOSOMES = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
		"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};

	/**
	 * Usage: java ExtendPhasingUsingRBP.java -g <path the genetic phasing vcf file> -p <path to the physical phasing vcf file>
	 * @param args -g <path the genetic phasing vcf file> -p <path to the physical phasing vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ExtendPhasingUsingRBP.java -g <path the genetic phasing vcf file> -p <path to the physical phasing vcf file>");
			System.exit(-1);
		} else {
			File geneticPhasingFile = null;
			File physicalPhasingFile = null;
			for (int i = 0; i <= 2; i += 2) {
				if (args[i].equals("-g")) {
					geneticPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					physicalPhasingFile = new File(args[i + 1]);
				}
			}
			try {
				extendPhasingUsingRBP(geneticPhasingFile, physicalPhasingFile);
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
	 * Tries to improve the result of a genetic phasing (eg: haplotyping, transmission phasing) 
	 * using the result of a physical phasing (ie read back phasing). 
	 * @param geneticPhasingFile partly phased vcf file generated using a genetic phasing method (haploscript or transmission)
	 * @param physicalPhasingFile partly phased vcf file generated using a physical phasing method (read backed phasing)
	 * @throws IOException
	 */
	private static void extendPhasingUsingRBP(File geneticPhasingFile, File physicalPhasingFile) throws IOException {
		// load genetic phasing vcf file		
		PhasedVectorList geneticVectorList = new PhasedVectorList();
		geneticVectorList.loadFromVCFFile(geneticPhasingFile);

		// load physical phasing vcf file
		PhasedVectorList physicalVectorList = new PhasedVectorList();
		physicalVectorList.loadFromVCFFile(physicalPhasingFile);

		mergeVectors(geneticVectorList, physicalVectorList);
	}


	/**
	 * Mergest the result of the genetic phasing and the physical phasing to try to get longer phasing 
	 * @param geneticVectorList
	 * @param physicalVectorList
	 */
	private static void mergeVectors(PhasedVectorList geneticVectorList, PhasedVectorList physicalVectorList) {
		// we analyze each chromosome defined in the CHROMOSOME constant
		for (String chromosome: CHROMOSOMES) {
			List<PhasedVector> geneticVectors = geneticVectorList.getPhasedVectorList(chromosome);
			List<PhasedVector> physicalVectors = physicalVectorList.getPhasedVectorList(chromosome);
			// we analyze each member of the family
			for (QuartetMember currentMember: QuartetMember.values()) {
				List<Integer> unphasedIndexes = new ArrayList<>();
				// if the genetic vectors and the physical vectors don't have the same size there is a pb with the files (or the program...)
				if ((geneticVectors != null) && (physicalVectors != null)) {
					if (geneticVectors.size() != physicalVectors.size()) {
						System.err.println("The genectic and physical lists don't have the same ammount of data");
					} else {
						// we create a list containing the index of all the genetic vectors not phased (that we can potentially phase)
						for (int i = 0; i < geneticVectors.size(); i++) {
							if (!geneticVectors.get(i).isPhased(currentMember)) {
								unphasedIndexes.add(i);
							}
						}
						int lastAnalyzedIndex = -1;
						// we analyze all the unphased variants one by one
						for (int currentUnphasedIndex: unphasedIndexes) {
							// we make sure that the current unphase variant has not been analyzed (
							if (currentUnphasedIndex > lastAnalyzedIndex) {
								if (physicalVectors.get(currentUnphasedIndex).isPhased(currentMember)) {
									int firstIndex = getFirstIndex(currentMember, currentUnphasedIndex, geneticVectors, physicalVectors);
									int lastIndex = getLastIndex(currentMember, currentUnphasedIndex, geneticVectors, physicalVectors);
									// if the first and last indexes are equals it means that we can't phase the genetic vector
									if (firstIndex != lastIndex) {
										List<PhasedVector> geneticVectorsToAnalyze = geneticVectors.subList(firstIndex, lastIndex);
										List<PhasedVector> physicalVectorsToAnalyze = physicalVectors.subList(firstIndex, lastIndex);
										PhasedVectorsInsert insert = new PhasedVectorsInsert(chromosome, currentUnphasedIndex, currentMember);
										insert.generateInsert(geneticVectorsToAnalyze, physicalVectorsToAnalyze);
										lastAnalyzedIndex = insert.getIndexLastPhasedVector();
									}
								}
							}
						}
					}
				}
			}
		}
	}


	/**
	 * @param member member of the quartet
	 * @param currentUnphasedIndex first unphased index
	 * @param geneticVectors list of genetic phasing vectors
	 * @param physicalVectors list of physical phasing vectors
	 * @return the index of the first mergeable vector of the genetic and physcal phasing
	 */
	private static int getFirstIndex(QuartetMember member, int currentUnphasedIndex, List<PhasedVector> geneticVectors, List<PhasedVector> physicalVectors) {
		boolean isPhased = true;
		// we go back as long as both the genetic and the physical vectors are phased and we don't reach the beginning of the lists
		while (isPhased && currentUnphasedIndex > 0) {
			isPhased = (physicalVectors.get(currentUnphasedIndex - 1).isPhased(member) && geneticVectors.get(currentUnphasedIndex - 1).isPhased(member));
			if (isPhased) {
				currentUnphasedIndex--;
			}
		}
		return currentUnphasedIndex;
	}


	/**
	 * @param member member of the quartet
	 * @param currentUnphasedIndex first unphased index
	 * @param geneticVectors list of genetic phasing vectors
	 * @param physicalVectors list of physical phasing vectors
	 * @return the index of the first mergeable vector of the genetic and physcal phasing
	 */	
	private static int getLastIndex(QuartetMember member, int currentUnphasedIndex, List<PhasedVector> geneticVectors, List<PhasedVector> physicalVectors) {
		boolean isPhased = true;
		boolean stillUnphased = true;
		while (isPhased && currentUnphasedIndex < physicalVectors.size() - 2) {			
			// still unphased is true as long as the genetic vectors are not phased.
			stillUnphased = stillUnphased && !geneticVectors.get(currentUnphasedIndex).isPhased(member);
			// isPhased is true if:
			// we are still in the unphased region of the genetic vectors and the physical vectors are phased  
			// we are back in a phased region and genetic and physical vectors are phased
			isPhased = (stillUnphased && physicalVectors.get(currentUnphasedIndex + 1).isPhased(member)) || 
					(!stillUnphased && physicalVectors.get(currentUnphasedIndex + 1).isPhased(member) && geneticVectors.get(currentUnphasedIndex + 1).isPhased(member));			
			if (isPhased) {
				currentUnphasedIndex++;
			}
		}
		return currentUnphasedIndex;
	}
}
