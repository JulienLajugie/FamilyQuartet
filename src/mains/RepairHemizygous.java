package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.QuartetMember;
import dataStructures.SegmentalDuplicationList;
import dataStructures.Variant;


/**
 * Generates a VCF file where the the hemizygous variants are repaired (0/. or 1/. instead of 0/0 or 1/1)
 * @author Julien Lajugie
 */
public class RepairHemizygous {


	/**
	 * Usage: java RepairHemizygous.java -v <path to the VCF file> -hp <path to the paternal hemizygous blocks file> 
	 * -hm <path to the maternal hemizygous blocks file> -hk1 <path to the kid1 hemizygous blocks file> -hk2 <path to the kid2 hemizygous blocks file>
	 * @param args -v <path to the VCF file> -hp <path to the paternal hemizygous blocks file> 
	 * -hm <path to the maternal hemizygous blocks file> -hk1 <path to the kid1 hemizygous blocks file> -hk2 <path to the kid2 hemizygous blocks file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java RepairHemizygous.java -v <path to the VCF file> -hp <path to the paternal hemizygous blocks file> " +
					"-hm <path to the maternal hemizygous blocks file> -hk1 <path to the kid1 hemizygous blocks file> -hk2 <path to the kid2 hemizygous blocks file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File paternalHemiBlockFile = null;
				File maternalHemiBlockFile = null;
				File kid1HemiBlockFile = null;
				File kid2HemiBlockFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					}
					if (args[i].equals("-hp")) {
						paternalHemiBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-hm")) {
						maternalHemiBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-hk1")) {
						kid1HemiBlockFile = new File(args[i + 1]);
					}
					if (args[i].equals("-hk2")) {
						kid2HemiBlockFile = new File(args[i + 1]);
					}
				}
				repairHemizygous(VCFFile, paternalHemiBlockFile, maternalHemiBlockFile, kid1HemiBlockFile, kid2HemiBlockFile);
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
		if (args.length != 10) {
			return false;
		}
		String[] mandatoryParameters = {"-v", "-hp", "-hm", "-hk1", "-hk2"};
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
	 * Generates a VCF file where the the hemizygous variants are repaired (0/. or 1/. instead of 0/0 or 1/1)
	 * @param VCFFile vcf file to repair
	 * @param paternalHemiBlockFile file with the paternal homzygous deletions
	 * @param maternalHemiBlockFile file with the maternal homzygous deletions
	 * @param kid1HemiBlockFile file with the kid1 homzygous deletions
	 * @param kid2HemiBlockFile file with the kid2 homzygous deletions
	 * @throws IOException
	 */
	private static void repairHemizygous(File VCFFile, File paternalHemiBlockFile, File maternalHemiBlockFile, File kid1HemiBlockFile, File kid2HemiBlockFile) throws IOException {
		SegmentalDuplicationList paternalHemiBlocks = new SegmentalDuplicationList();
		paternalHemiBlocks.loadBedOrBgr(paternalHemiBlockFile);
		SegmentalDuplicationList maternalHemiBlocks = new SegmentalDuplicationList();
		maternalHemiBlocks.loadBedOrBgr(maternalHemiBlockFile);
		SegmentalDuplicationList kid1HemiBlocks = new SegmentalDuplicationList();
		kid1HemiBlocks.loadBedOrBgr(kid1HemiBlockFile);
		SegmentalDuplicationList kid2HemiBlocks = new SegmentalDuplicationList();
		kid2HemiBlocks.loadBedOrBgr(kid2HemiBlockFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					//System.out.println(line);
				} else {
					String correctedLine = line;
					try {
						Variant currentVariant = new Variant(line);
						if (paternalHemiBlocks.getBlock(currentVariant.getChromosome(), currentVariant.getPosition()) != null) { 
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.FATHER);							
						}
						if (maternalHemiBlocks.getBlock(currentVariant.getChromosome(), currentVariant.getPosition()) != null) { 
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.MOTHER);							
						}
						if (kid1HemiBlocks.getBlock(currentVariant.getChromosome(), currentVariant.getPosition()) != null) { 
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.KID1);							
						}
						if (kid2HemiBlocks.getBlock(currentVariant.getChromosome(), currentVariant.getPosition()) != null) { 
							correctedLine = getCorrectedLine(correctedLine, QuartetMember.KID2);							
						}
						if (correctedLine == null) {
							//currentVariant.printVariantBgrFormat();
							// if we weren't able to repair the line we just print the line with no modification
							correctedLine = line;
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
	 * @return the corrected line
	 */
	private static String getCorrectedLine(String lineToCorrect, QuartetMember member) {
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
		if (splitLine[indexGenotypeInfo].charAt(0) == splitLine[indexGenotypeInfo].charAt(2)) {
			splitLine[indexGenotypeInfo] = "." + splitLine[indexGenotypeInfo].substring(1);
		} else {
			return null;
		}
		for (int i = 1; i < splitLine.length; i++) {
			correctedLine += '\t' + splitLine[i];
		}
		return correctedLine;
	}
}
