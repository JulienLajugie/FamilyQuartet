package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.Variant;
import exceptions.FilteredVCFLineException;
import exceptions.VCFException;


/**
 * Adds a PLF flags to the info field of the VCF file when on the VCF line rejected by the PL filter
 * @author Julien Lajugie
 */
public class MarkPLFilteredVariants {

	private static final String PLF_INFO_HEADER = "##INFO=<ID=PLF,Number=0,Type=Flag,Description=\"Variant rejected by the " +
			"PL Filter: min(PLs) < "+ Variant.INDIVIDUALS_PL_MIN_VALUE + "\">";	// header for the RDF flag
	private static final int	INFO_FIELD_INDEX = 7; // index of the info field in the vcf file

	/**
	 * Usage: java MarkPLFilteredVariants -v <path to the VCF file>
	 * @param args -v <path to the VCF file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java MarkPLFilteredVariants.java -v <path to the VCF file>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-v")) {
						VCFFile = new File(args[i + 1]);
					}
				}
				markPLFilteredVariants(VCFFile);
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
		// case with no -v parameter
		if (!args[0].equals("-v")) {
			return false;
		}
		return true;
	}


	/**
	 * Adds a PLF flags to the info field of the VCF file when on the VCF line rejected by the PL filter
	 * @param VCFFile
	 * @throws IOException
	 */
	private static void markPLFilteredVariants(File VCFFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFile));
			String line = null;
			boolean infoStarted = false;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) == '#') {
					if (line.startsWith("##INFO")) {
						infoStarted = true;
					} else if (infoStarted) {
						System.out.println(PLF_INFO_HEADER);
						infoStarted= false;
					}
					System.out.println(line);
				} else {
					try {
						new Variant(line);
					} catch (FilteredVCFLineException e) {
						// if the variant got rejected by the PL filter we mark it
						if (e.getFilterName().equalsIgnoreCase("PL")) {
							line = markLineAs(line, "PLF");
						}
					} catch (VCFException e) {
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
	 * @param line a vcf line
	 * @param flag a flag to add to the info field
	 * @return a copy of the vcf line where the specified flag was added to the info field 
	 */
	private static String markLineAs(String line, String flag) {
		String[] splitLine = line.split("\t");
		splitLine[INFO_FIELD_INDEX] += ";" + flag; 
		String newLine = splitLine[0];
		for (int i = 1; i < splitLine.length; i++) {
			newLine += "\t" + splitLine[i];
		}
		return newLine;
	}
}
