package dataStructures;

import java.util.ArrayList;
import java.util.List;

/**
 * This class represents an insertion of readba
 * @author Julien Lajugie
 */
public class PhasedVectorsInsert {
	
	private List<PhasedVector> 	insertPhasedVectors; 			// list of consecutive vectors phased with physical phasing to insert in genetic phasing
	private final String		chromosome;						// chromosome of the insert
	private final int 			indexFirstInsertVector;			// index of the element in the genetic and physical list corresponding to the first vector of the insert 
	private int 				indexLastInsertVector;			// index of the element in the genetic and physical list corresponding to the last vector of the insert
	private final QuartetMember	quartetMember;					// quartet member of the insert
	private int					goodVectorBeforeInsertCount;	// number of physical vectors having a phasing compatible with the genetic vectors before the insert 
	private int					badVectorBeforeInsertCount;		// number of physical vectors having a phasing incompatible with the genetic vectors before the insert
	private int					goodVectorAfterInsertCount;		// number of physical vectors having a phasing compatible with the genetic vectors after the insert 
	private int					badVectorAfterInsertCount;		// number of physical vectors having a phasing incompatible with the genetic vectors after the insert
	private boolean				isGeneticVectorInvertionNeeded;	// true if all the genetic vectors after the junction needs to be inverted (cross over happened in kid1)
	
	
	/**
	 * Creates an instance of {@link PhasedVectorsInsert} 
	 * @param chromosome chromosome of the insert
	 * @param indexFirstInsertVector index of the element in the genetic and physical list of the first vector of the insert
	 * @param quartetMember quartet member of the insert
	 */
	public PhasedVectorsInsert(String chromosome, int indexFirstInsertVector, QuartetMember quartetMember) {
		this.chromosome = chromosome;
		this.indexFirstInsertVector = indexFirstInsertVector;
		this.quartetMember = quartetMember;
	}
	
	
	/**
	 * Generates the insert as well as the the statistics about the insert (error and supporting vectors)
	 * @param geneticList list of vectors partially phased using a genetic method (haplotyping or transmission phasing) 
	 * @param physicalList list of vectors  phased using a physical method (read backed phasing)
	 */
	public void generateInsert(List<PhasedVector> geneticList, List<PhasedVector> physicalList) {
		// search the index of the first unphased genetic vector (can be zero)
		// this index correspond to the first element of the index
		int indexFirstUnphased = 0;
		while ((indexFirstUnphased < geneticList.size()) && (geneticList.get(indexFirstUnphased).isPhased(quartetMember))) {
			indexFirstUnphased++;
		}
		// search the index of the last unphased genetic vector (can be the last index)
		// this index correspond to the last element of the insert
		int indexLastUnphased = indexFirstUnphased;
		while ((indexLastUnphased < geneticList.size()) && (!geneticList.get(indexLastUnphased).isPhased(quartetMember))) {
			indexLastUnphased++;
		}
		if (geneticList.get(indexLastUnphased).isPhased(quartetMember)) {
			indexLastUnphased--;
		}
		this.indexLastInsertVector = this.indexFirstInsertVector + (indexLastUnphased - indexFirstUnphased);
		this.insertPhasedVectors = physicalList.subList(indexFirstUnphased, indexLastUnphased);
		// compute how many vectors are identical and how many vectors are not identical before the insert
		for (int i = 0; i < indexFirstUnphased; i++) {
			if (geneticList.get(i).getGenotype(quartetMember).equals(physicalList.get(i).getGenotype(quartetMember))) {
				goodVectorBeforeInsertCount++;
			} else {
				badVectorBeforeInsertCount++;
			}
		}
		// compute how many vectors are identical and how many vectors are not identical after the insert
		for (int i = indexLastUnphased + 1; i < geneticList.size(); i++) {
			if (geneticList.get(i).getGenotype(quartetMember).equals(physicalList.get(i).getGenotype(quartetMember))) {
				goodVectorAfterInsertCount++;
			} else {
				badVectorAfterInsertCount++;
			}
		}
		analyzeInser();		
	}
	
	
	
	private void analyzeInser() {
		/*
		 * TODO create the method that tells if:
		 *  - the insert is valid or not
		 *  - the insert needs to be inverted
		 *  - the rest of the genetic phasing needs to be inverted
		 */
		
	}
	
	
	/**
	 * Inverts the phasing of the insert vector for the studied member
	 */
	private void invertInsert() {
		for (PhasedVector currentVector: insertPhasedVectors) {
			currentVector.invert(quartetMember);
		}
	}


	/**
	 * @return the list of consecutive vectors phased with physical phasing to insert in genetic phasing
	 */
	public final List<PhasedVector> getInsertPhasedVectors() {
		return insertPhasedVectors;
	}
	
	
	/**
	 * @return the chromosome of the insert
	 */
	public final String getChromosome() {
		return chromosome;
	}
	
	
	/**
	 * @return the index of the element in the genetic and physical list corresponding to the first vector of the insert
	 */
	public final int getIndexFirstPhasedVector() {
		return indexFirstInsertVector;
	}
	
	
	/**
	 * @return the index of the element in the genetic and physical list corresponding to the last vector of the insert
	 */
	public final int getIndexLastPhasedVector() {
		return indexLastInsertVector;
	}
	
	
	/**
	 * @return the quartet member of the insert
	 */
	public final QuartetMember getQuartetMember() {
		return quartetMember;
	}
	
	
	/**
	 * @return the number of physical vectors having a phasing compatible with the genetic vectors before the insert
	 */
	public final int getGoodVectorBeforeInsertCount() {
		return goodVectorBeforeInsertCount;
	}
	
	
	/**
	 * @return the number of physical vectors having a phasing incompatible with the genetic vectors before the insert
	 */
	public final int getBadVectorBeforeInsertCount() {
		return badVectorBeforeInsertCount;
	}
	

	/**
	 * @return the number of physical vectors having a phasing compatible with the genetic vectors after the insert
	 */
	public final int getGoodVectorAfterInsertCount() {
		return goodVectorAfterInsertCount;
	}
	
	
	/**
	 * @return the number of physical vectors having a phasing incompatible with the genetic vectors after the insert
	 */
	public final int getBadVectorAfterInsertCount() {
		return badVectorAfterInsertCount;
	}
	
		
	/**
	 * @return true if all the genetic vectors after the junction needs to be inverted (cross over happened in kid1)
	 */
	public final boolean isGeneticVectorInvertionNeeded() {
		return isGeneticVectorInvertionNeeded;
	}
}
