package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.CrossTriosInheritanceState;
import dataStructures.InheritanceStateBlockList;
import dataStructures.InheritanceStateBlockListFactory;
import dataStructures.Variant;


/**
 * Generates a VCF file where the contamination from the mother to the dad has been corrected 
 * (causing false ab+aa;ab/ab that are in fact aa/bb;ab/ab)
 * @author Julien Lajugie
 */
public class RepairContaminatedVariants {


	/**
	 * Usage: java RepairContaminatedVariants -v <path to the VCF file> -b <block file (optional)>
	 * @param args -v <path to the VCF file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java RepairContaminatedVariants.java -v <path to the VCF file> -b <block file (optional)>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File blockFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					}
					if (args[i].equals("-b")) {
						blockFile = new File(args[i + 1]);
					}	
				}
				repairContaminatedVariants(VCFFile, blockFile);
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
		if ((args.length != 2) && (args.length != 4)) {
			return false;
		}
		if (args[0].equals("-v")) {
			return true;
		} else {
			return false;
		}
	}


	/**
	 * Generates a VCF file where the contamination from the mother to the dad has been corrected 
	 * (causing false ab+aa;ab/ab that are in fact aa/bb;ab/ab)
	 * @param VCFFile VCF files with the variants of the family quartet
	 * @param blockFile file containing the inheritance state block. can be null
	 * @throws IOException if the VCF file is not valid
	 */
	private static void repairContaminatedVariants(File VCFFile, File blockFile) throws IOException {
		InheritanceStateBlockList<CrossTriosInheritanceState> blockList = null;
		if (blockFile != null) {
			blockList = InheritanceStateBlockListFactory.createFromCrossTriosBgrFile(blockFile);
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					System.out.println(line);
				} else {
					String correctedLine = line;
					try {
						Variant currentVariant = new Variant(line);
						// if the block file was specified we just want to correct the MIE
						if ((blockList == null) || 
								((blockList.getBlock(currentVariant) != null) 
										&& (blockList.getBlock(currentVariant).isSCE(currentVariant)))) {
							if ((currentVariant.getGenotypePattern().equals("ab+aa;ab/ab") 
									|| currentVariant.getGenotypePattern().startsWith("ab/ab"))
									&& !currentVariant.getGenotypePattern().endsWith("bb")) {
								String fatherRepairedGenotype = currentVariant.getContaminationCorrectedGenotype(line);
								if (fatherRepairedGenotype != null) {
									correctedLine = getCorrectedVcfLine(line, fatherRepairedGenotype);
									// if the correction created a MIE we discard it
									Variant newVariant = new Variant(correctedLine);
									if (newVariant.isMIE()) {
										correctedLine = line;
									}
									// if the correction created a SCE we discard it
									if ((blockList != null) && (blockList.getBlock(currentVariant) != null) && blockList.getBlock(currentVariant).isSCE(currentVariant)) {
										correctedLine = line;
									}
								}
							}
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
	 * @param line input vcf line before correction
	 * @param fatherRepairedGenotype corrected paternal genotype
	 * @return the vcf line with the corrected paternal genotype
	 */
	private static String getCorrectedVcfLine(String line, String fatherRepairedGenotype) {
		String[] splitLine = line.split("\t");
		// retrieve the paternal genotype info field from the vcf
		String paternalGenotype = splitLine[9].trim();
		String paternalGenotypeInfoFields[] = splitLine[9].trim().split(":");
		// replace the paternal genotype info field
		paternalGenotype = fatherRepairedGenotype;
		// add paternal info fields until PL field
		paternalGenotype += ':' + paternalGenotypeInfoFields[1] + ':' + paternalGenotypeInfoFields[2] + ':' + paternalGenotypeInfoFields[3] + ':';
		// replace PL field with arbitrary score
		if (fatherRepairedGenotype.charAt(0) == '0') {
			paternalGenotype += "0,99,99";
		} else {
			paternalGenotype += "99,99,0";
		}
		// add remaining genotype info fields
		for (int i = 5; i < paternalGenotypeInfoFields.length; i++) {
			paternalGenotype += ':' + paternalGenotypeInfoFields[i];
		}
		// recreate the vcfLine
		String phasedVcfLine = "";
		for (int i = 0; i < 9; i++) {
			phasedVcfLine += splitLine[i] + "\t";
		}
		phasedVcfLine += paternalGenotype + "\t";
		phasedVcfLine += splitLine[10].trim() +'\t' + splitLine[11].trim() + '\t' + splitLine[12].trim();
		return phasedVcfLine;
	}
}
