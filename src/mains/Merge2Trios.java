package mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


/**
 * Merges the phasing of 2 trios phased by transmission in a phased quartet file
 * @author Julien Lajugie
 */
public class Merge2Trios {

	/**
	 * Usage: java Merge2Trios.java -t1 <path to the 1st trio VCF> -t2 <path to the 2nd trio VCF>
	 * @param args -t1 <path to the 1st trio VCF> -t2 <path to the 2nd trio VCF>
	 */
	public static void main(String[] args) {
		// exit the program if the input parameters are not correct
		if (!areParametersValid(args)) { 
			System.out.println("Usage: java Merge2Trios.java -t1 <path to the 1st trio VCF> -t2 <path to the 2nd trio VCF>");
			System.exit(-1);
		} else {
			File trio1File = null;
			File trio2File = null;
			for (int i = 0; i <= 2; i += 2) {
				if (args[i].equals("-t1")) {
					trio1File = new File(args[i + 1]);
				}
				if (args[i].equals("-t2")) {
					trio2File = new File(args[i + 1]);
				}
			}
			try {
				merge2Trios(trio1File, trio2File);
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
		// case with no -t1 parameter
		if (!args[0].equals("-t1") && !args[2].equals("-t1")) {
			return false;
		}
		// case with no -t2 parameter
		if (!args[0].equals("-t2") && !args[2].equals("-t2")) {
			return false;
		}
		return true;
	}


	/**
	 * Merges the phasing of 2 trios phased by transmission in a phased quartet file
	 * @param trio1File 1st trio vcf file
	 * @param trio2File 2nd trio vcf file
	 * @throws IOException
	 */
	private static void merge2Trios(File trio1File, File trio2File) throws IOException {
		BufferedReader readerTrio1 = null;
		BufferedReader readerTrio2 = null;
		try {
			readerTrio1 = new BufferedReader(new FileReader(trio1File));
			readerTrio2 = new BufferedReader(new FileReader(trio2File));
			String lineTrio1 = null;
			String lineTrio2 = null;
			// loop until eof
			while (((lineTrio1 = readerTrio1.readLine()) != null) && ((lineTrio2 = readerTrio2.readLine()) != null)){
				// a line starting with a # is a comment line
				if (lineTrio1.charAt(0) == '#') {
					System.out.println(lineTrio1);
				} else {
					String[] splitLineTrio1 = lineTrio1.split("\t");
					String[] splitLineTrio2 = lineTrio2.split("\t");
					if (!splitLineTrio1[0].equals(splitLineTrio2[0])) {
						throw new IOException("Files cannot be merged");
					} else {
						String outputLine = "";
						int currentFied = 0;
						while (currentFied < splitLineTrio1.length - 1) {
							outputLine += splitLineTrio1[currentFied] + '\t';
							currentFied++;							
						}
						outputLine += splitLineTrio2[currentFied];
						System.out.println(outputLine);
					}					
				}
			}
		} finally {
			if (readerTrio1 != null) {
				readerTrio1.close();
			}
			if (readerTrio2 != null) {
				readerTrio2.close();
			}
		}
	}
}
