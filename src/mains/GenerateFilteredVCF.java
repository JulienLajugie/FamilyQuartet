package mains;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.InheritanceState;
import dataStructures.Variant;
import exceptions.InvalidVCFLineException;


/**
 * Generates a filtered VCF file.  The filters must be defined in the Variant class.
 * @author Julien Lajugie
 */
public class GenerateFilteredVCF {
	
	
	/**
	 * Usage: java GenerateFilteredVCF -f <path to the file>
	 * @param args -f <path to the file>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if	((args.length != 2) ||
				(!args[0].equals("-f"))) {
			System.out.println("Usage: java GenerateFilteredVCF -f <path to the file>");
			System.exit(-1);
		} else {
			try {
				generateFilteredVCF(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * Generates a filtered VCF file.  The filters must be defined in the Variant class.
	 * @param VCFFilePath VCF files with the variants
	 * @throws IOException if the VCF file is not valid
	 */
	private static void generateFilteredVCF(String VCFFilePath) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(VCFFilePath));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// a line starting with a # is a comment line
				if (line.charAt(0) != '#') {
					try {
						Variant currentVariant = new Variant(line);
						// we don't process variants with more than one alternative allele or indels
						if ((currentVariant.getAlternatievAllele().length() != 1) || (currentVariant.getReferenceAllele().length() != 1)) {
							throw new InvalidVCFLineException();
						}
						// we don't process variants that are not informative or a mandelien inheritance state
						if ((currentVariant.getInheritanceStates()[0] == InheritanceState.NOT_INFORMATIVE) &&
								(currentVariant.getInheritanceStates()[0] == InheritanceState.MIE)) {
							throw new InvalidVCFLineException();
						}
						System.out.println(line);
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
}
