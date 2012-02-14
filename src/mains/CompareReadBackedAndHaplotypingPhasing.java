package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.PhasedGenotypesSeries;
import dataStructures.PhasedVector;
import dataStructures.PhasedVectorList;


/**
 * Compares the result between a read back phasing and a haplotyping phasing
 * @author Julien Lajugie
 */
public class CompareReadBackedAndHaplotypingPhasing {

	/**
	 * Usage: java CompareReadBackedAndHaplotypingPhasing -h <path to haplotyping phased vector file> -v <path to the vcf file>
	 * @param args -h <path to haplotyping phased vector file> -v <path to the vcf file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java CompareReadBackedAndHaplotypingPhasing -h <path to haplotyping phased vector file> -v <path to the vcf file>");
			System.exit(-1);
		} else {
			File haplotypingFile = null;
			File vcfFile = null;
			for (int i = 0; i <= 2; i += 2) {
				if (args[i].equals("-h")) {
					haplotypingFile = new File(args[i + 1]);
				}
				if (args[i].equals("-v")) {
					vcfFile = new File(args[i + 1]);
				}
			}
			try {
				compareReadBackedAndHaplotypingPhasing(haplotypingFile, vcfFile);
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
		// case with no -h parameter
		if (!args[0].equals("-h") && !args[2].equals("-h")) {
			return false;
		}
		// case with no -v parameter
		if (!args[0].equals("-v") && !args[2].equals("-v")) {
			return false;
		}
		return true;
	}


	/**
	 * Compares the result between a read back phasing and a haplotyping phasing
	 * @param haplotypingFile haplotyping file
	 * @param vcfFile vcf file
	 * @throws IOException 
	 */
	private static void compareReadBackedAndHaplotypingPhasing(File haplotypingFile, File vcfFile) throws IOException {
		// load haplotyping file		
		PhasedVectorList haplotypingVectorList = new PhasedVectorList();
		//haplotypingVectorList.loadFromHaplotypingFile(haplotypingFile);
		haplotypingVectorList.loadFromVCFFile(haplotypingFile);
		// load vcf file
		PhasedVectorList vcfVectorList = new PhasedVectorList();
		vcfVectorList.loadFromVCFFile(vcfFile);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(haplotypingFile));
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
					PhasedVector haplotypingPhasedVector = haplotypingVectorList.getPhasedVector(chromosome, position);
					PhasedVector vcfPhasedVector = vcfVectorList.getPhasedVector(chromosome, position);
					// father
					if ((haplotypingPhasedVector == null) || (vcfPhasedVector == null)) {
						//System.out.println(chromosome + '\t' + position + '\t' + haplotypingPhasedVector + '\t' + vcfPhasedVector);
					} else {
						int result = paternalSeries.addGenotypes(haplotypingPhasedVector.getFatherGenotype(), vcfPhasedVector.getFatherGenotype());
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += paternalSeries.getCompatibleGenotypes();
							incompatibleGenotypesCount += paternalSeries.getIncompatibleGenotypes();
							if (paternalSeries.getIncompatibleGenotypes() > 0) {
								paternalSeries.print();
							}
							paternalSeries.reset();
						}
						// mother
						result = maternalSeries.addGenotypes(haplotypingPhasedVector.getMotherGenotype(), vcfPhasedVector.getMotherGenotype());
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += maternalSeries.getCompatibleGenotypes();
							incompatibleGenotypesCount += maternalSeries.getIncompatibleGenotypes();
							maternalSeries.reset();
						}
						// kid 1
						result = kid1Series.addGenotypes(haplotypingPhasedVector.getKid1Genotype(), vcfPhasedVector.getKid1Genotype());
						if (result == PhasedGenotypesSeries.SERIES_FINISHED) {
							compatibleGenotypesCount += kid1Series.getCompatibleGenotypes();
							incompatibleGenotypesCount += kid1Series.getIncompatibleGenotypes();
							kid1Series.reset();
						}
						// kid 2
						result = kid2Series.addGenotypes(haplotypingPhasedVector.getKid2Genotype(), vcfPhasedVector.getKid2Genotype());
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
