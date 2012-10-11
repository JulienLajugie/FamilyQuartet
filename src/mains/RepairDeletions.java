package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.SegmentalDuplicationList;


/**
 * Generates a VCF file by adding the deletions from the result in bed files
 * @author Julien Lajugie
 */
public class RepairDeletions {

	/**
	 * VCF info header for the hemizygous error filter flag
	 */
	private static final String HEF_INFO_HEADER = "##INFO=<ID=HEF,Number=0,Type=Flag,Description=\"Hemizygous Error Filter\">";
	private static final int	INFO_FIELD_INDEX = 7; // index of the info field in the vcf file

	/**
	 * Usage: java RepairDeletions.java -v <path to the VCF file> -dp <path to the paternal deleted blocks file> 
	 * -dm <path to the maternal deleted blocks file> -dk1 <path to the kid1 deleted blocks file> -dk2 <path to the kid2 deleted blocks file> -n <number of deleted allele>
	 * @param args -v <path to the VCF file> -dp <path to the paternal deleted blocks file> 
	 * -dm <path to the maternal deleted blocks file> -dk1 <path to the kid1 deleted blocks file> -dk2 <path to the kid2 deleted blocks file> -n <number of deleted allele>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			showUsage();
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File paternalDeletedBlockFile = null;
				File maternalDeletedBlockFile = null;
				File kid1DeletedBlockFile = null;
				File kid2DeletedBlockFile = null;
				int alleleDeletedCount = 0;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					}
					if (args[i].equals("-dp")) {
						paternalDeletedBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-dm")) {
						maternalDeletedBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-dk1")) {
						kid1DeletedBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-dk2")) {
						kid2DeletedBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-n")) {
						alleleDeletedCount = Integer.parseInt(args[i + 1]);
					}
				}
				if ((alleleDeletedCount != 1) && (alleleDeletedCount != 2)) {
					System.out.println("Incorrect allele count");
					showUsage();
					System.exit(-1);
				}
				repairDeletions(VCFFile, paternalDeletedBlockFile, maternalDeletedBlockFile, kid1DeletedBlockFile, kid2DeletedBlockFile, alleleDeletedCount);
			} catch (Exception e) {
				e.printStackTrace();
				showUsage();
			}
		}
	}


	/**
	 * Shows the usage
	 */
	private static void showUsage() {
		System.out.println("Usage: java RepairDeletions.java -v <path to the VCF file> -dp <path to the paternal deleted blocks file>" +
				" -dm <path to the maternal deleted blocks file> -dk1 <path to the kid1 deleted blocks file> -dk2 <path to the kid2 deleted blocks file> -n <number of deleted allele>");
	}


	/**
	 * @param args parameters from the main function
	 * @return true if the parameters are valid
	 */
	private static boolean areParametersValid(String[] args) {
		if (args == null) {
			return false;
		}
		if (args.length != 12) {
			return false;
		}
		String[] mandatoryParameters = {"-v", "-dp", "-dm", "-dk1", "-dk2", "-n"};
		for (String currentMandatoryParameter: mandatoryParameters) {
			boolean found = false;
			int i = 0;
			while ((i < args.length) && !found) {
				found = args[i].equals(currentMandatoryParameter); 
				i += 2;
			}
			if (!found) {
				return false;
			}
		}
		return true;
	}


	/**
	 * Generates a VCF file by adding the deletions from the result in bed files
	 * @param VCFFile input VCF file
	 * @param paternalDeletedBlockFile file with the paternal deletion blocks
	 * @param maternalDeletedBlockFile file with the maternal deletion blocks
	 * @param kid1DeletedBlockFile file with the kid1 deletion blocks
	 * @param kid2DeletedBlockFile file with the kid2 deletion blocks
	 * @param alleleDeletedCount number of allele that got deleted
	 * @throws IOException
	 */
	private static void repairDeletions(File VCFFile, File paternalDeletedBlockFile, File maternalDeletedBlockFile, File kid1DeletedBlockFile, File kid2DeletedBlockFile, int alleleDeletedCount) throws IOException {
		SegmentalDuplicationList paternalHemiBlocks = new SegmentalDuplicationList();
		paternalHemiBlocks.loadBedOrBgr(paternalDeletedBlockFile);
		SegmentalDuplicationList maternalHemiBlocks = new SegmentalDuplicationList();
		maternalHemiBlocks.loadBedOrBgr(maternalDeletedBlockFile);
		SegmentalDuplicationList kid1HemiBlocks = new SegmentalDuplicationList();
		kid1HemiBlocks.loadBedOrBgr(kid1DeletedBlockFile);
		SegmentalDuplicationList kid2HemiBlocks = new SegmentalDuplicationList();
		kid2HemiBlocks.loadBedOrBgr(kid2DeletedBlockFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			boolean infoStarted = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					// add the hemizygous error filter flag to the vcf file in the case of hemizygous corrections
					if (alleleDeletedCount == 1) {
						if (line.startsWith("##INFO")) {
							infoStarted = true;
						} else if (infoStarted) {
							System.out.println(HEF_INFO_HEADER);
							infoStarted= false;
						}
					}
					System.out.println(line);
				} else {
					String[] splitLine = line.split("\t");
					String correctedLine = line;
					try {
						String chromosome = splitLine[0].trim();
						int position = Integer.parseInt(splitLine[1].trim());
						if (paternalHemiBlocks.getBlock(chromosome, position) != null) {
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.FATHER, alleleDeletedCount);
						}
						if (maternalHemiBlocks.getBlock(chromosome, position) != null) {
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.MOTHER, alleleDeletedCount);
						}
						if (kid1HemiBlocks.getBlock(chromosome, position) != null) {
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.KID1, alleleDeletedCount);
						}
						if (kid2HemiBlocks.getBlock(chromosome, position) != null) {
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.KID2, alleleDeletedCount);
						}
					} catch (Exception e) {
						// do nothing
					} finally {
						System.out.println(correctedLine);
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
	 * @param lineToCorrect the line that needs to be corrected
	 * @param member the member that needs to be corrected
	 * @param alleleDeletedCount number of allele that got deleted
	 * @return the corrected line
	 */
	private static String getCorrectedLine(String lineToCorrect, QuartetMember member,int alleleDeletedCount) {
		if (lineToCorrect == null) { 
			return null;
		}
		int indexGenotypeInfo = 0;
		String[] splitLine = lineToCorrect.split("\t");
		String correctedLine = splitLine[0];
		switch (member) {
		case FATHER:
			indexGenotypeInfo = 9;
			break;
		case MOTHER:
			indexGenotypeInfo = 10;
			break;
		case KID1:
			indexGenotypeInfo = 11;
			break;
		case KID2:
			indexGenotypeInfo = 12;
			break;
		}
		if (alleleDeletedCount == 2) {
			splitLine[indexGenotypeInfo] = ".|." + splitLine[indexGenotypeInfo].substring(3);
		} else {
			assert alleleDeletedCount == 1 : "incorrect deleted allele count";

			if (splitLine[indexGenotypeInfo].charAt(0) == splitLine[indexGenotypeInfo].charAt(2)) {
				splitLine[indexGenotypeInfo] = "." + splitLine[indexGenotypeInfo].substring(1);
			} else {
				// if it's not already an hemizygous and the genotype is not homozygous we flag the line as HEF
				if (splitLine[indexGenotypeInfo].charAt(0) != '.' && splitLine[indexGenotypeInfo].charAt(2) != '.') {
					lineToCorrect = markLineAs(lineToCorrect, "HEF");
					return lineToCorrect;
				}
			}
		}
		for (int i = 1; i < splitLine.length; i++) {
			correctedLine += '\t' + splitLine[i];
		}
		return correctedLine;
	}


	/**
	 * @param line a vcf line
	 * @param flag a flag to add to the info field
	 * @return a copy of the vcf line where the specified flag was added to the info field if it wasn't already present
	 */
	private static String markLineAs(String line, String flag) {
		String[] splitLine = line.split("\t");
		if (splitLine[INFO_FIELD_INDEX].contains(flag)) {
			return line;
		}
		splitLine[INFO_FIELD_INDEX] += ";" + flag;
		String newLine = splitLine[0];
		for (int i = 1; i < splitLine.length; i++) {
			newLine += "\t" + splitLine[i];
		}
		return newLine;
	}
}
