package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


/**
 * Repairs the output phased vector file of the software Haploscript by completing the incomplete vectors 
 * @author Julien Lajugie
 */
public class RepairHaploscriptOutput {


	/**
	 * Usage: java RepairHaploscriptOutput -f <path to haploscript phased vector file>
	 * @param args -f <path to haploscript phased vector file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java RepairHaploscriptOutput -f <path to haploscript phased vector file>");
			System.exit(-1);
		} else {
			File phasedVectorFile = null;
			phasedVectorFile = new File(args[1]);
			try {
				repairHaploscriptOutput(phasedVectorFile);
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
		if (args.length != 2) {
			return false;
		}
		// case with no -f parameter
		if (!args[0].equals("-f")) {
			return false;
		}
		return true;
	}


	/**
	 * Repairs the output phased vector file of the software Haploscript by completing the incomplete vectors 
	 * @param phasedVectorFile phased vector file from Haplotyping 
	 * @throws IOException
	 */
	private static void repairHaploscriptOutput(File phasedVectorFile) throws IOException {

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(phasedVectorFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					System.out.println(line);
				} else {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					String position = splitLine[1].trim();
					String binaryState = splitLine[2].trim();
					String unphasedVector = splitLine[3].trim();
					String phasedVector = splitLine[4].trim();
					String correctedPhasedVector = "";
					if (binaryState.contains(".")) {
						correctedPhasedVector = phasedVector;
					} else {
						correctedPhasedVector = correctVector(unphasedVector, phasedVector);
					}
					String newVcfLine = chromosome + '\t' + position + '\t' + binaryState + '\t' + unphasedVector + '\t' + correctedPhasedVector;
					System.out.println(newVcfLine);
				}
			}	
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Corrects a phased vector from the Haploscript software
	 * @param unphasedVector
	 * @param phasedVector
	 * @return
	 */
	private static String correctVector(String unphasedVector, String phasedVector) {
		if ((phasedVector.equals("........")) || (phasedVector.length() != 8))  {
			return phasedVector;
		}
		char paternal1stHaplotype = unphasedVector.charAt(0);
		char paternal2ndHaplotype = unphasedVector.charAt(1);
		char maternal1stHaplotype = unphasedVector.charAt(2);
		char maternal2ndHaplotype = unphasedVector.charAt(3);
		char kid1PaternalAllele = phasedVector.charAt(4);
		char kid1MaternalAllele = phasedVector.charAt(5);		
		String correctedVector = "";
		correctedVector +=  kid1PaternalAllele;
		correctedVector += (paternal1stHaplotype != kid1PaternalAllele ? paternal1stHaplotype : paternal2ndHaplotype);
		correctedVector += kid1MaternalAllele;
		correctedVector += (maternal1stHaplotype != kid1MaternalAllele ? maternal1stHaplotype : maternal2ndHaplotype);
		correctedVector += phasedVector.substring(4);
		return correctedVector;
	}
}
