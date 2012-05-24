package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.PhasedGenotypesSeries;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;


/**
 * Compares the result between a physical phasing and a genetic phasing
 * @author Julien Lajugie
 */
public class ComparePhysicalAndGeneticPhasing {

	/**
	 * Usage: java ComparePhysicalAndGeneticPhasing.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>
	 * @param args -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java ComparePhysicalAndGeneticPhasing.java -g <path to genetic phasing vcf file> -p <path to physical phasing vcf file>");
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
				comparePhysicalAndGeneticPhasing(geneticPhasingFile, physicalPhasingFile);
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
		if (!args[0].equals("-") && !args[2].equals("-p")) {
			return false;
		}
		return true;
	}


	/**
	 * Compares the result between a physical phasing and a genetic phasing
	 * @param geneticPhasingFile file with the result of the genetic phasing
	 * @param physicalPhasingFile file with the result of the physical phasing
	 * @throws IOException 
	 */
	private static void comparePhysicalAndGeneticPhasing(File geneticPhasingFile, File physicalPhasingFile) throws IOException {
		// load genetic phasing file		
		PhasedVectorList geneticVectorList = new PhasedVectorList();
		geneticVectorList.loadFromVCFFile(geneticPhasingFile);

		// load physical phasing file
		PhasedVectorList physicalVectorList = new PhasedVectorList();
		physicalVectorList.loadFromVCFFile(physicalPhasingFile);
		physicalVectorList.convertFromPhysicalToGeneticVectorList();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(geneticPhasingFile));
			String line = null;
			int compatibleGenotypesCount = 0; // count of variants with a compatible read back and haplotyping phasing
			int incompatibleGenotypesCount = 0; // count of variants with a incompatible read back and haplotyping phasing
			PhasedGenotypesSeries paternalSeries = new PhasedGenotypesSeries();
			PhasedGenotypesSeries maternalSeries = new PhasedGenotypesSeries();
			PhasedGenotypesSeries kid1Series = new PhasedGenotypesSeries();
			PhasedGenotypesSeries kid2Series = new PhasedGenotypesSeries();
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					PhasedVector geneticPhasedVector = geneticVectorList.getPhasedVector(chromosome, position);
					PhasedVector physicalPhasedVector = physicalVectorList.getPhasedVector(chromosome, position);
					// father
					if ((geneticPhasedVector != null) && (physicalPhasedVector != null)) {
						int result = paternalSeries.addGenotypes(geneticPhasedVector.getFatherGenotype(), physicalPhasedVector.getFatherGenotype(), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += paternalSeries.getCompatibleGenotypes();
							incompatibleGenotypesCount += paternalSeries.getIncompatibleGenotypes();
							paternalSeries.reset();
						}
						// mother
						result = maternalSeries.addGenotypes(geneticPhasedVector.getMotherGenotype(), physicalPhasedVector.getMotherGenotype(), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += maternalSeries.getCompatibleGenotypes();
							incompatibleGenotypesCount += maternalSeries.getIncompatibleGenotypes();
							maternalSeries.reset();
						}
						// kid 1
						result = kid1Series.addGenotypes(geneticPhasedVector.getKid1Genotype(), physicalPhasedVector.getKid1Genotype(), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += kid1Series.getCompatibleGenotypes();
							incompatibleGenotypesCount += kid1Series.getIncompatibleGenotypes();
							kid1Series.reset();
						}
						// kid 2
						result = kid2Series.addGenotypes(geneticPhasedVector.getKid2Genotype(), physicalPhasedVector.getKid2Genotype(), chromosome, position);
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += kid2Series.getCompatibleGenotypes();
							incompatibleGenotypesCount += kid2Series.getIncompatibleGenotypes();							
							kid2Series.reset();
						}
					}
				}
			}
			double ratioCompatible = compatibleGenotypesCount / (double) (compatibleGenotypesCount + incompatibleGenotypesCount) * 100d; 
			System.out.println("Compatible vector count: " + compatibleGenotypesCount
					+ ", incompatible vector count: " + incompatibleGenotypesCount
					+ ", ratio: " + ratioCompatible);
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}
}
