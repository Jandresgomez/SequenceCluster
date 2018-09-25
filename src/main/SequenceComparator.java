package main;

import java.util.Comparator;
import java.util.stream.IntStream;

import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;

public class SequenceComparator implements Comparator<DNAShortKmer> {

	@Override
	public int compare(DNAShortKmer s1, DNAShortKmer s2) {
		int k = 0;
		int l = Math.min(s1.length(),s2.length());
		
		for(int i = 0; i < l; i++) {
			if(s1.charAt(i) - s2.charAt(i) != 0) k++;
		}
		
		return k;
	}

}
