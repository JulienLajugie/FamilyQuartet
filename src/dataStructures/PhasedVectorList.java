package dataStructures;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import exceptions.InvalidVCFFieldException;
import exceptions.VCFException;


/**
 * This class represent a list of phased vector loadable from a Haploscript file
 * @author Julien Lajugie
 */
public class PhasedVectorList {

	private final Map<String, List<PhasedVector>> phasedVectorMap; // list of phased vector organised by chromosome


	/**
	 * Creates an instance of {@link PhasedVectorList}
	 */
	public PhasedVectorList() {
		phasedVectorMap = new HashMap<String, List<PhasedVector>>();
	}


	/**
	 * Load a phased vector file
	 * @param phasedVectorFile
	 * @throws IOException
	 */
	public void loadFromHaplotypingFile(File phasedVectorFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(phasedVectorFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't care about the comment lines 
				if (line.trim().charAt(0) != '#') {
					String[] splitLine = line.split("\t");
					String chromosome = splitLine[0].trim();
					int position = Integer.parseInt(splitLine[1].trim());
					String unphasedVector = splitLine[3].trim();
					String phasedVector = splitLine[4].trim();
					if ((!phasedVector.equals("SCE")) && (!phasedVector.equals("MIE"))) {
						PhasedVector vectorToAdd = new PhasedVector(position, unphasedVector, phasedVector);
						//System.out.println(chromosome + '\t' + position + '\t' + unphasedVector +'\t' + phasedVector + '\t' + vectorToAdd.getFatherGenotype() + '\t' + vectorToAdd.getMotherGenotype() + '\t' + vectorToAdd.getKid1Genotype() + '\t' + vectorToAdd.getKid2Genotype());
						// if the list doesn't contain the chromosome we add it
						if (!phasedVectorMap.containsKey(chromosome)) {
							List<PhasedVector> listToAdd = new ArrayList<PhasedVector>();
							listToAdd.add(vectorToAdd);
							phasedVectorMap.put(chromosome, listToAdd);						
						} else {
							phasedVectorMap.get(chromosome).add(vectorToAdd);
						}
					}
				}
			}
			// we sort the list in position order
			sortLists();
		} finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * Load a vcf File
	 * @param vcfFile a vcf file
	 * @throws IOException
	 */
	public void loadFromVCFFile(File vcfFile) throws IOException {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(vcfFile));
			String line = null;
			// loop until eof
			while ((line = reader.readLine()) != null) {
				// we don't care about the comment lines 
				if (line.trim().charAt(0) != '#') {
					try {
						Variant variant = new Variant(line);
						int position = variant.getPosition();
						String chromosome = variant.getChromosome();
						// generate the phased vector
						String phasedVector = "";
						phasedVector += alleleToPhasedVector(variant.getFatherAlleles()[0], variant.isFatherPhased());
						phasedVector += alleleToPhasedVector(variant.getFatherAlleles()[1], variant.isFatherPhased());
						phasedVector += alleleToPhasedVector(variant.getMotherAlleles()[0], variant.isMotherPhased());
						phasedVector += alleleToPhasedVector(variant.getMotherAlleles()[1], variant.isMotherPhased());
						phasedVector += alleleToPhasedVector(variant.getKid1Alleles()[0], variant.isKid1Phased());
						phasedVector += alleleToPhasedVector(variant.getKid1Alleles()[1], variant.isKid1Phased());
						phasedVector += alleleToPhasedVector(variant.getKid2Alleles()[0], variant.isKid2Phased());
						phasedVector += alleleToPhasedVector(variant.getKid2Alleles()[1], variant.isKid2Phased());
						// generate the unphased vector
						String unphasedVector = "";
						unphasedVector += alleleToUnphasedVector(variant.getFatherAlleles()[0]);
						unphasedVector += alleleToUnphasedVector(variant.getFatherAlleles()[1]);
						unphasedVector += alleleToUnphasedVector(variant.getMotherAlleles()[0]);
						unphasedVector += alleleToUnphasedVector(variant.getMotherAlleles()[1]);
						unphasedVector += alleleToUnphasedVector(variant.getKid1Alleles()[0]);
						unphasedVector += alleleToUnphasedVector(variant.getKid1Alleles()[1]);
						unphasedVector += alleleToUnphasedVector(variant.getKid2Alleles()[0]);
						unphasedVector += alleleToUnphasedVector(variant.getKid2Alleles()[1]);
						PhasedVector vectorToAdd = new PhasedVector(position, unphasedVector, phasedVector);
						//System.out.println(chromosome + '\t' + position + '\t' + unphasedVector +'\t' + phasedVector + '\t' + vectorToAdd.getFatherGenotype() + '\t' + vectorToAdd.getMotherGenotype() + '\t' + vectorToAdd.getKid1Genotype() + '\t' + vectorToAdd.getKid2Genotype());
						// if the list doesn't contain the chromosome we add it
						if (!phasedVectorMap.containsKey(chromosome)) {
							List<PhasedVector> listToAdd = new ArrayList<PhasedVector>();
							listToAdd.add(vectorToAdd);
							phasedVectorMap.put(chromosome, listToAdd);
						} else {
							phasedVectorMap.get(chromosome).add(vectorToAdd);
						}
					} catch (VCFException e) {
						if (e instanceof InvalidVCFFieldException) {
							e.printStackTrace();
							System.err.println(((InvalidVCFFieldException) e).getVCFFieldName() + " - " + ((InvalidVCFFieldException) e).getVCFFieldValue());
						}
					}				
				}
			}
			// we sort the list in position order
			sortLists();
		}  finally {
			if (reader != null) {
				reader.close();
			}
		}
	}


	/**
	 * @param allele an {@link AlleleType}
	 * @param isPhased true if the sample is phased
	 * @return a phased vector item for the specified allele
	 */
	private char alleleToPhasedVector(AlleleType allele, boolean isPhased) {
		if ((allele == null) || (!isPhased)) {
			return '.';
		}
		return (allele == AlleleType.REFERENCE_ALLELE ? 'a' : 'b');
	}


	/**
	 * @param allele an {@link AlleleType}
	 * @return an unphased vector item for the specified allele
	 */	
	private char alleleToUnphasedVector(AlleleType allele) {
		if (allele == null) {
			return '.';
		}
		return (allele == AlleleType.REFERENCE_ALLELE ? 'a' : 'b');
	}



	/**
	 * For each chromosome of the map, this method sorts the {@link PhasedVector} by position
	 */
	public void sortLists() {
		for (List<PhasedVector> currentList: phasedVectorMap.values()) {
			Collections.sort(currentList);
		}
	}


	/**
	 * @param chromosome a chromosome
	 * @param position a position
	 * @return the {@link PhasedVector} on the specified chromosome at the specified position
	 */
	public PhasedVector getPhasedVector(String chromosome, int position) {
		if (!phasedVectorMap.containsKey(chromosome)) {
			return null;
		}
		List<PhasedVector> vectorList = phasedVectorMap.get(chromosome);
		return findPhasedVector(vectorList, position);		
	}


	/**
	 * Binary search of {@link PhasedVector} sorted by position.
	 * @param list list of phased vectors
	 * @param position position of the phased vector to find
	 * @return the {@link PhasedVector} at the specified position.  Null it there is none
	 */
	private PhasedVector findPhasedVector(List<PhasedVector> list, int position) {
		int start = 0;
		int stop = list.size() - 1;
		while (start < stop) {
			int mid = (start + stop) / 2;  // Compute mid point.
			if (position < list.get(mid).getPosition()) {
				stop = mid;  // repeat search in bottom half.
			} else if (position > list.get(mid).getPosition()) {
				start = mid + 1;  // Repeat search in top half.
			} else {
				return list.get(mid); // Found it. return position
			}
		}
		return null; // Failed to find key
	}
}
