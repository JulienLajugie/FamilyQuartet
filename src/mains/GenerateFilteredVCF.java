package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import dataStructures.SegmentalDuplicationList;
import dataStructures.Variant;
import exceptions.VCFException;


/**
 * Generates a filtered VCF file.  The filters must be defined in the Variant class.
 * @author Julien Lajugie
 */
public class GenerateFilteredVCF {


	/**
	 * Usage: java GenerateFilteredVCF -f <path to the file> -s <segmental duplication file (optional)>
	 * @param args -f <path to the file> -s <segmental duplication file (optional)>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java GenerateFilteredVCF -f <path to the file> -s <segmental duplication file (optional)>");
			System.exit(-1);
		} else {
			try {
				File VCFFile = null;
				File segDupFile = null;
				for (int i = 0; i < args.length; i += 2) {
					if (args[i].equals("-f")) {
						VCFFile = new File(args[i + 1]);
					} else if (args[i].equals("-s")) {			
						segDupFile = new File(args[i + 1]);
					}					
				}
				generateFilteredVCF(VCFFile, segDupFile);
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
		if (args[0].equals("-f") || args[2].equals("-f")) {
			return true;
		} else {
			return false;
		}
	}


	/**
	 * Generates a filtered VCF file.  The filters must be defined in the Variant class.
	 * @param VCFFile VCF files with the variants
	 * @param segDupFile bed or bgr file containing the segmental duplication. Variants in these regions will be excluded.  Can be null
	 * @throws IOException if the VCF file is not valid
	 */
	private static void generateFilteredVCF(File VCFFile, File segDupFile) throws IOException {
		SegmentalDuplicationList segDupList = null;
		if (segDupFile != null) {
			segDupList = new SegmentalDuplicationList(segDupFile);
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
					try {
						Variant currentVariant = new Variant(line);
						if ((segDupList != null) && !segDupList.isInSegmentalDuplication(currentVariant)) {
							if (!currentVariant.isIndel()) {
								System.out.println(line);
							}
						}
					} catch (VCFException e) {
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
