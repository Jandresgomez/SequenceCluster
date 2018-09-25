package main;

import java.util.HashMap;

import ngsep.sequences.DNAShortKmer;

public class ClusterTableManager {
	private int newIndex = 0;
	private int[][] table;
	private HashMap<DNAShortKmer, Integer> index;
	private static int KMER_LENGTH = 31;
	private static final char[] ALPH = {'A', 'G', 'C', 'T'};
	
	public ClusterTableManager(int size) {
		table = new int[4][size*KMER_LENGTH];
		index = new HashMap<>(size/2);
	}
	
	public void add(DNAShortKmer kmer) {
		boolean hit = false;
		for(int i = 0; i < kmer.length(); i++) {
			for(int j = 0; j < ALPH.length; j++) {
				Integer k = index.get(new DNAShortKmer(kmer.subSequence(0, i) + "" + ALPH[j] + kmer.subSequence(i+1, kmer.length())));
				if(k != null) {
					hit = true;
					append(kmer, k);
				}
			}
		}
		
		if(!hit) {
			newCluster(kmer);
		}
	}
	
	public void newCluster(DNAShortKmer kmer) {
		for(int i = 0; i < kmer.length(); i++) {
			for(int j = 0; j < ALPH.length; j++) {
				if(kmer.charAt(i) == ALPH[j]) { table[j][newIndex*KMER_LENGTH + i]++; break; }
			}
		}
		
		index.put(kmer, newIndex++);
	}
	
	public DNAShortKmer getRepresentative(Integer k) {
		char[] structure = new char[KMER_LENGTH];
		
		for(int i = 0; i < structure.length; i++) {
			int max = 0;
			for(int j = 0; j < ALPH.length; j++) {
				if(max <= table[j][k*KMER_LENGTH + i])
					structure[i] = ALPH[j];
			}
		}
		
		String seq = "";
		for(char c : structure) seq += c;
		return new DNAShortKmer(seq);
	}
	
	public void append(DNAShortKmer kmer, Integer k) {
		DNAShortKmer oldKmer = getRepresentative(k);
		
		for(int i = 0; i < kmer.length(); i++) {
			for(int j = 0; j < ALPH.length; j++) {
				if(kmer.charAt(i) == ALPH[j]) { table[j][k*KMER_LENGTH + i]++; break; }
			}
		}
		
		DNAShortKmer newKmer = getRepresentative(k);
		index.remove(oldKmer);
		index.put(newKmer, k);
	}
	
	public static void main(String[] args) {
		DNAShortKmer kmer = new DNAShortKmer("ACGT"); 
		for(int i = 0; i < kmer.length(); i++) {
			for(int j = 0; j < ALPH.length; j++) {
				//System.out.println(kmer.subSequence(0, i) + "===" + kmer.subSequence(i+1, kmer.length()));
				System.out.println(new DNAShortKmer(kmer.subSequence(0, i) + "" + ALPH[j] + kmer.subSequence(i+1, kmer.length())));
			}
		}
	}
}
