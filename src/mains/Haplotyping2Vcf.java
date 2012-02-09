package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.InheritanceStateBlock;
import dataStructures.InheritanceStateBlockList;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;

public class Haplotyping2Vcf {

	private static final String PHASED_SET_FORMAT_HEADER = 
			"##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase Set\">";
	private static final String PHASED_SET_FORMAT_FIELD = "PS";

	
	/**
	 * Usage: java Haplotyping2Vcf -v <path to the VCF file> -b <path to the block file> -p <path to the phased vector file>
	 * @param args -v <path to the VCF file> -b <path to the block file> -p <path to the phased vector file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java Haplotyping2Vcf -v <path to the VCF file> -b <path to the block file> -p <path to the phased vector file>");
			System.exit(-1);
		} else {
			File VCFFile = null, blockFile = null, phasedVectorFile = null;
			for (int i = 0; i <= 4; i += 2) {
				if (args[i].equals("-v")) {
					VCFFile = new File(args[i + 1]);
				}
				if (args[i].equals("-b")) {
					blockFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					phasedVectorFile = new File(args[i + 1]);
				}
			}
			try {
				haplotyping2Vcf(VCFFile, blockFile, phasedVectorFile);
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
		if (!args[0].equals("-p") && !args[2].equals("-p") && !args[4].equals("-p")) {
			return false;
		}
		return true;

	}


	/**
	 * Generates a phased VCF file filtered for ISCA
	 * @param VCFFile input VCF file
	 * @param blockFile ISCA file with the inheritance state blocks
	 * @param phasedVectorFile phased vector file from Haplotyping 
	 * @throws IOException
	 */
	private static void haplotyping2Vcf(File VCFFile, File blockFile, File phasedVectorFile) throws IOException {
		// load the inheritance block list
		InheritanceStateBlockList blockList = new InheritanceStateBlockList();
		blockList.loadFromISCAFile(blockFile);
		// load the phased vector file
		PhasedVectorList vectorList = new PhasedVectorList();
		vectorList.loadFromFile(phasedVectorFile);

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			boolean formatHeaderSet = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					if ((line.contains("##FORMAT")) && !formatHeaderSet) {
						// print the header for the 2 new fields
						System.out.println(PHASED_SET_FORMAT_HEADER);
						formatHeaderSet = true;
						//System.out.println(BLOCK_INFO_FORMAT_HEADER);
					}
					System.out.println(line);
				} else {
					Variant variant;
					try {
						String[] splitLine = line.split("\t");
						String chromosome = splitLine[0].trim();
						int position = Integer.parseInt(splitLine[1].trim());
						InheritanceStateBlock block = blockList.getBlock(chromosome, position);
						variant = new Variant(line);
						String newVcfLine = "";
						if (variant.getGenotypePattern().equals("aa/aa;aa/aa")) {
							// fully homozygous vectors are not phased by haploscripting
							newVcfLine = substituteVcfLine(splitLine, block, variant);
						} else {
							PhasedVector phasedVector = vectorList.getPhasedVector(chromosome, position);
							newVcfLine = substituteVcfLine(splitLine, block, phasedVector);
						}
						System.out.println(newVcfLine);
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


	/**
	 * @param splitLine vcf line of the current variant split in fields delimited by tabulation 
	 * @param block block containing the current variant.  null if none
	 * @param phasedVector vector phased by haploscripting
	 * @return the phased vcf line
	 */
	private static String substituteVcfLine(String[] splitLine, InheritanceStateBlock block, PhasedVector phasedVector) {
		String formatField = splitLine[8].trim();
		String fatherGenotype = splitLine[9].trim();
		String motherGenotype = splitLine[10].trim();
		String kid1Genotype = splitLine[11].trim();
		String kid2Genotype = splitLine[12].trim();
		if (phasedVector != null) {			
			fatherGenotype = substituteGenotype(fatherGenotype, phasedVector.getFatherGenotype());
			motherGenotype = substituteGenotype(motherGenotype, phasedVector.getMotherGenotype());
			kid1Genotype = substituteGenotype(kid1Genotype, phasedVector.getKid1Genotype());
			kid2Genotype = substituteGenotype(kid2Genotype, phasedVector.getKid2Genotype());
			if (block != null) {
				String phaseSet = Integer.toString(block.getStartPosition());
				formatField += ":" + PHASED_SET_FORMAT_FIELD;
				fatherGenotype += ":" + phaseSet;
				motherGenotype += ":" + phaseSet;
				kid1Genotype += ":" + phaseSet;
				kid2Genotype += ":" + phaseSet;
			}
		}
		String phasedVcfLine = "";
		for (int i = 0; i < 8; i++) {
			phasedVcfLine += splitLine[i] + "\t";
		}
		phasedVcfLine += formatField + "\t";
		phasedVcfLine += fatherGenotype + "\t";
		phasedVcfLine += motherGenotype + "\t";
		phasedVcfLine += kid1Genotype + "\t";
		phasedVcfLine += kid2Genotype;
		return phasedVcfLine;
	}


	/**
	 * @param vcfGenotype vcf genotype field
	 * @param phasedGenotype phased genotype
	 * @return a phased genotype field of a vcf file
	 */
	private static String substituteGenotype(String vcfGenotype, String phasedGenotype) {
		String[] splitVcfGenotype = vcfGenotype.split(":");
		String resultGenotype = phasedGenotype;
		for (int i = 1; i < splitVcfGenotype.length; i++) {
			resultGenotype += ":" + splitVcfGenotype[i];
		}
		return resultGenotype;
	}


	/**
	 * @param splitLine vcf line of the current variant split in fields delimited by tabulation 
	 * @param block block containing the current variant.  null if none
	 * @param Variant current variant
	 * @return the phased vcf line of a homozygous variant
	 */
	private static String substituteVcfLine(String[] splitLine, InheritanceStateBlock block, Variant variant) {
		String formatField = splitLine[8].trim();
		String fatherGenotype = splitLine[9].trim();
		String motherGenotype = splitLine[10].trim();
		String kid1Genotype = splitLine[11].trim();
		String kid2Genotype = splitLine[12].trim();
		if (variant != null) {			
			fatherGenotype = fatherGenotype.replace('/', '|');
			motherGenotype = motherGenotype.replace('/', '|');
			kid1Genotype = kid1Genotype.replace('/', '|');
			kid2Genotype = kid2Genotype.replace('/', '|');
			if (block != null) {
				String phaseSet = Integer.toString(block.getStartPosition());
				formatField += ":" + PHASED_SET_FORMAT_FIELD;
				fatherGenotype += ":" + phaseSet;
				motherGenotype += ":" + phaseSet;
				kid1Genotype += ":" + phaseSet;
				kid2Genotype += ":" + phaseSet;
			}
		}
		String phasedVcfLine = "";
		for (int i = 0; i < 8; i++) {
			phasedVcfLine += splitLine[i] + "\t";
		}
		phasedVcfLine += formatField + "\t";
		phasedVcfLine += fatherGenotype + "\t";
		phasedVcfLine += motherGenotype + "\t";
		phasedVcfLine += kid1Genotype + "\t";
		phasedVcfLine += kid2Genotype;
		return phasedVcfLine;
	}
}
