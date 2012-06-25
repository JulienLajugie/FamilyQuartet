package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.PhasedGenotypesSeries;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;
import dataStructures.QuartetMember;


/**
 * Generates a bgr file where for the specified member a position gets a score of 1 if the physical 
 * and the genetic phasing are compatible and a score of -1 if they are not  
 * @author Julien Lajugie
 */
public class GeneratePhysicalVsGeneticPhasingBgr {

	/**
	 * Usage: java GeneratePhysicalVsGeneticPhasingBgr.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -m <quartet member>
	 * @param args -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -m <quartet member>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GeneratePhysicalVsGeneticPhasingBgr.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -m <quartet member>");
			System.exit(-1);
		} else {
			File geneticPhasingFile = null;
			File physicalPhasingFile = null;
			QuartetMember member = null;
			for (int i = 0; i < args.length; i += 2) {
				if (args[i].equals("-g")) {
					geneticPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-p")) {
					physicalPhasingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-m")) {
					member = QuartetMember.valueOf(args[i + 1]);
					if (member == null) {
						System.out.println("Usage: java GeneratePhysicalVsGeneticPhasingBgr.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file> -m <quartet member>");
						System.exit(-1);						
					}
				}
			}
			try {
				generatePhysicalVsGeneticPhasingBgr(geneticPhasingFile, physicalPhasingFile, member);
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
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m") && !args[4].equals("-m")) {
			return false;
		}
		return true;
	}


	/**
	 * Generates a bgr file where for the specified member a position gets a score of 1 if the physical 
	 * and the genetic phasing are compatible and a score of -1 if they are not  
	 * @param geneticPhasingFile file with the result of the genetic phasing
	 * @param physicalPhasingFile file with the result of the physical phasing
	 * @param member a member of the family quartet
	 * @throws IOException 
	 */
	private static void generatePhysicalVsGeneticPhasingBgr(File geneticPhasingFile, File physicalPhasingFile, QuartetMember member) throws IOException {
		// load genetic phasing file		
		PhasedVectorList geneticVectorList = new PhasedVectorList();
		geneticVectorList.loadFromVCFFile(geneticPhasingFile);

		// load physical phasing file
		PhasedVectorList physicalVectorList = new PhasedVectorList();
		physicalVectorList.loadFromVCFFile(physicalPhasingFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(physicalPhasingFile));
			String line = null;
			PhasedGenotypesSeries phasedSeries = new PhasedGenotypesSeries();
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
						int result = phasedSeries.addGeneticPhysicalGenotypes(geneticPhasedVector.getGenotype(member), physicalPhasedVector.getGenotype(member), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							phasedSeries.printResultPhasingBgr();
							phasedSeries.reset();
						}
					}
				}
			}
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
