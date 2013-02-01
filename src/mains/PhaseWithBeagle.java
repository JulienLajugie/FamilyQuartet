package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.SegmentalDuplication;
import dataStructures.SegmentalDuplicationList;
import dataStructures.TrioInheritanceState;
import dataStructures.Variant;

/**
 * Uses the result from the beagle phasing to phase parental blocks and name the allele in the VCF file
 * @author Julien Lajugie
 */
public class PhaseWithBeagle {

	private static String AA_HEADER = "##FORMAT=<ID=AA,Number=2,Type=String,Description=\"Ancestral Alleles\">"; // header for ancestral allele


	/**
	 * Usage: java PhaseWithBeagle -v <path to the VCF file> -b <path to the block file> -p <path to the beagle paternal block file> -m <path to the beagle maternal block file>
	 * @param args -v <path to the VCF file> -b <path to the block file> -p <path to the beagle paternal block file> -m <path to the beagle maternal block file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) {
			System.out.println("Usage: java PhaseWithBeagle -v <path to the VCF file> -b <path to the block file> -p <path to the beagle paternal block file> -m <path to the beagle maternal block file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File blockFile = null;
				File paternalBlockFile = null;
				File maternalBlockFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					} else if (args[i].equals("-b")) {
						blockFile = new File(args[i + 1]);
					} else if (args[i].equals("-p")) {
						paternalBlockFile = new File(args[i + 1]);
					} else if (args[i].equals("-m")) {
						maternalBlockFile = new File(args[i + 1]);
					}
				}
				phaseWithBeagle(VCFFile, blockFile, paternalBlockFile, maternalBlockFile);
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
		if (args.length != 8) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v") && !args[4].equals("-v") && !args[6].equals("-v")) {
			return false;
		}
		// case with no -b parameter
		if (!args[0].equals("-b") && !args[2].equals("-b") && !args[4].equals("-b") && !args[6].equals("-b")) {
			return false;
		}
		// case with no -p parameter
		if (!args[0].equals("-p") && !args[2].equals("-p") && !args[4].equals("-p") && !args[6].equals("-p")) {
			return false;
		}
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m") && !args[4].equals("-m") && !args[6].equals("-m")) {
			return false;
		}
		return true;
	}


	/**
	 * Uses the result from the beagle phasing to phase parental blocks and name the allele in the VCF file
	 * @param vcfFile
	 * @param blockFile
	 * @param paternalBlockFile
	 * @param maternalBlockFile
	 */
	private static void phaseWithBeagle(File vcfFile, File blockFile, File paternalBlockFile, File maternalBlockFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> ISBlockList;
		ISBlockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);

		SegmentalDuplicationList paternalBeagleBlocks = new SegmentalDuplicationList();
		paternalBeagleBlocks.loadBedOrBgrWithScore(paternalBlockFile);

		SegmentalDuplicationList maternalBeagleBlocks = new SegmentalDuplicationList();
		maternalBeagleBlocks.loadBedOrBgrWithScore(maternalBlockFile);

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			String currentChromo = null;
			String previousChromo = null;
			SegmentalDuplication previousPaternalBeagleBlock = null;
			SegmentalDuplication previousMaternalBeagleBlock = null;
			String kid1PA = null; // kid1 paternal allele
			String kid1MA = null; // kid1 maternal allele
			boolean formatStarted = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					if (!line.startsWith("##FORMAT=<ID=PS")) { // we want to remove the phase set field
						if (line.startsWith("##FORMAT")) {
							formatStarted = true;
						} else if (formatStarted) {
							System.out.println(AA_HEADER);
							formatStarted= false;
						}
						System.out.println(line);
					}
				} else {
					try {
						Variant variant = new Variant(line);
						InheritanceStateBlock<CrossTriosInheritanceState> ISBlock = ISBlockList.getBlock(variant);
						if ((ISBlock != null) && !variant.isMIE() && !variant.isSCE(ISBlock.getBlockState())) {
							currentChromo = variant.getChromosome();
							if (!currentChromo.equals(previousChromo)) {
								// reset the previsous beagle blocks when we change chromosome
								previousPaternalBeagleBlock = null;
								previousMaternalBeagleBlock = null;
								// a chromosome starts with the kid1 having the 1st parental alleles
								kid1PA = "P1";
								kid1MA = "M1";
							}
							// check if the kid1 paternal allele changed
							SegmentalDuplication currentPaternalBeagleBlock = paternalBeagleBlocks.getBlock(currentChromo, variant.getPosition());
							if ((currentPaternalBeagleBlock != null) && !currentPaternalBeagleBlock.equals(previousPaternalBeagleBlock)) {
								if ((previousPaternalBeagleBlock != null) && (currentPaternalBeagleBlock.getScore() != previousPaternalBeagleBlock.getScore())) {
									kid1PA = invertParentalAllele(kid1PA);
								}
								previousPaternalBeagleBlock = currentPaternalBeagleBlock;
							}
							// check if the kid1 maternal allele changed
							SegmentalDuplication currentMaternalBeagleBlock = maternalBeagleBlocks.getBlock(currentChromo, variant.getPosition());
							if ((currentMaternalBeagleBlock != null) && !currentMaternalBeagleBlock.equals(previousMaternalBeagleBlock)) {
								if ((previousMaternalBeagleBlock != null) && (currentMaternalBeagleBlock.getScore() != previousMaternalBeagleBlock.getScore())) {
									kid1MA = invertParentalAllele(kid1MA);
								}
								previousMaternalBeagleBlock = currentMaternalBeagleBlock;
							}
							// compute the kid2 ancestral alleles
							String kid2PA = computeKid2ParentalAllele(kid1PA, ISBlock.getBlockState().getPaternalTrioState());
							String kid2MA = computeKid2ParentalAllele(kid1MA, ISBlock.getBlockState().getMaternalTrioState());
							// check if the dad genotype needs to be changed
							boolean invertDad = paternalBeagleBlocks.getBlock(variant.getChromosome(), variant.getPosition()).getScore() == -1;
							// check if the mom genotype needs to be changed
							boolean invertMom = maternalBeagleBlocks.getBlock(variant.getChromosome(), variant.getPosition()).getScore() == -1;
							line = addAncestralAllele(line, invertDad, invertMom, kid1PA, kid1MA, kid2PA, kid2MA);
							previousChromo = currentChromo;
						}
					} catch (Exception e) {
						//e.printStackTrace();
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
	 * @param kid1Parental
	 * @param parentalState
	 * @return the kid2 parental allele depending on kid1 and the inheritance state block
	 */
	private static String computeKid2ParentalAllele(String kid1Parental, TrioInheritanceState parentalState) {
		switch (parentalState) {
		case IDENTICAL:
			// in this case kid 2 as the same parental allele as kid1
			return kid1Parental;
		case UNKNOWN:
			// return P? or M?
			return kid1Parental.charAt(0) + "?";
		case NON_IDENTICAL:
			// in this case kid 2 as the opposite parental allele as kid1
			return invertParentalAllele(kid1Parental);
		default:
			// shouldn't happen
			return null;
		}
	}


	/**
	 * @param parentalAllele a parental allele
	 * @return the opposite parental allele
	 */
	private static String invertParentalAllele(String parentalAllele) {
		String invertedParentalAllele = "1";
		if (parentalAllele.charAt(1) == '1') {
			invertedParentalAllele = "2";
		}
		return parentalAllele.charAt(0) + invertedParentalAllele;
	}

	/**
	 * @param line vcf line
	 * @param isBlock block of the current variant
	 * @param variant variant of the vcf line
	 * @return the vcf line with the parental blocks phased and name the allele in the VCF file
	 */
	private static String addAncestralAllele(String line, boolean invertDad, boolean invertMom, String kid1PA, String kid1MA, String kid2PA, String kid2MA) {
		String[] splitLine = line.split("\t");
		// we remove the phase set from the GT info fields if it is present
		if (splitLine[8].contains("PS")) {
			for (int i = 8; i < splitLine.length; i++) {
				splitLine[i] = removePhaseSetFromGTInfoField(splitLine[i]);
			}
		} else {
			for (int i = 8; i < splitLine.length; i++) {
				splitLine[i] += ":";
			}
		}
		splitLine[8] += "AA";
		// we invert the dad genotype if needed and add the ancestral alleles field
		splitLine[9] += "P1,P2";
		if (invertDad) {
			splitLine[9] = invertGenotype(splitLine[9]);
		}
		// we invert the mom genotype if needed and add the ancestral alleles field
		splitLine[10] += "M1,M2";
		if (invertMom) {
			splitLine[10] = invertGenotype(splitLine[10]);
		}
		// add the kid1 ancestral alleles field
		splitLine[11] += kid1MA +"," + kid1PA;
		// add the kid2 ancestral alleles field
		splitLine[12] += kid2MA +"," + kid2PA;
		String newLine = splitLine[0];
		for (int i = 1; i < splitLine.length; i++) {
			newLine += "\t" + splitLine[i];
		}
		return newLine;
	}


	/**
	 * @param genotypeField
	 * @return the specified genotype info field with the genotype inverted
	 */
	private static String invertGenotype(String genotypeField) {
		char oldAllele1 = genotypeField.charAt(0);
		char phasing = genotypeField.charAt(1);
		char oldAllele2 = genotypeField.charAt(2);
		String invertedGT = Character.toString(oldAllele2) + Character.toString(phasing) + Character.toString(oldAllele1) + genotypeField.substring(3);
		return invertedGT;
	}


	/**
	 * @param gtInfoField
	 * @return the specified genotype info field with the phase set removed
	 */
	private static String removePhaseSetFromGTInfoField(String gtInfoField) {
		String[] splitInfo = gtInfoField.split(":");
		String newInfoField = "";
		for (int i = 0; i < (splitInfo.length -1); i++) {
			newInfoField += splitInfo[i] + ":";
		}
		return newInfoField;
	}
}
