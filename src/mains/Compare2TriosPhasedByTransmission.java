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
 * Compares the result between two trios (Father / Mother / Kid1 & Father / Mother / Kid2)
 * for a specified founder to find the inheritance state of each variant 
 * (founder identical / not identical)
 * @author Julien Lajugie
 */
public class Compare2TriosPhasedByTransmission {

	/**
	 * Usage: java Compare2TriosPhasedByTransmission.java -t1 <path to first trio vcf file> -t2 <path to second trio vcf file> -m <FATHER or MOTHER>
	 * @param args -t1 <path to first trio vcf file> -t2 <path to second trio vcf file> -m <FATHER or MOTHER>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java Compare2TriosPhasedByTransmission.java -t1 <path to first trio vcf file> -t2 <path to second trio vcf file> -m <FATHER or MOTHER>");
			System.exit(-1);
		} else {
			File vcfTrio1 = null;
			File vcfTrio2 = null;
			QuartetMember member = null;
			for (int i = 0; i <= 4; i += 2) {
				if (args[i].equals("-t1")) {
					vcfTrio1 = new File(args[i + 1]);
				}
				if (args[i].equals("-t2")) {
					vcfTrio2 = new File(args[i + 1]);
				}
				if (args[i].equals("-m")) {
					member = QuartetMember.valueOf(args[i + 1]);
					if ((member == null) || ((member != QuartetMember.FATHER) && (member != QuartetMember.MOTHER))) {
						System.out.println("Usage: java Compare2TriosPhasedByTransmission.java -t1 <path to first trio vcf file> -t2 <path to second trio vcf file> -m <FATHER or MOTHER>");
						System.exit(-1);						
					}
				}
			}
			try {
				compare2TriosPhasedByTransmission(vcfTrio1, vcfTrio2, member);
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
		// case with no -t1 parameter
		if (!args[0].equals("-t1") && !args[2].equals("-t1") && !args[4].equals("-t1")) {
			return false;
		}
		// case with no -t2 parameter
		if (!args[0].equals("-t2") && !args[2].equals("-t2") && !args[4].equals("-t2")) {
			return false;
		}
		// case with no -m parameter
		if (!args[0].equals("-m") && !args[2].equals("-m") && !args[4].equals("-m")) {
			return false;
		}
		return true;
	}



	/**
	 * Compares the result between two trios (Father / Mother / Kid1 & Father / Mother / Kid2)
	 * for a specified founder to find the inheritance state of each variant 
	 * (founder identical / not identical)
	 * @param vcfTrio1 vcf with the transmission phasing of the 1st trio 
	 * @param vcfTrio2 vcf with the transmission phasing of the 2nd trio
	 * @param quartetMember a founder {@link QuartetMember}
	 * @throws IOException
	 */
	private static void compare2TriosPhasedByTransmission(File vcfTrio1, File vcfTrio2, QuartetMember quartetMember) throws IOException {
		// load trio 1		
		PhasedVectorList trio1VectorList = new PhasedVectorList();
		trio1VectorList.loadFromVCFFile(vcfTrio1);
		// load trio 2
		PhasedVectorList trio2VectorList = new PhasedVectorList();
		trio2VectorList.loadFromVCFFile(vcfTrio2);		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfTrio1));
			String line = null;
			PhasedGenotypesSeries phasedVectorSeries = new PhasedGenotypesSeries();
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					PhasedVector trio1PhasedVector = trio1VectorList.getPhasedVector(chromosome, position);
					PhasedVector trio2PhasedVector = trio2VectorList.getPhasedVector(chromosome, position);
					if ((trio1PhasedVector != null) && (trio2PhasedVector != null)) {
						int result = phasedVectorSeries.add2GeneticGenotypes(trio1PhasedVector.getGenotype(quartetMember), trio2PhasedVector.getGenotype(quartetMember), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							phasedVectorSeries.printResultPhasingBgr();
							phasedVectorSeries.reset();
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
